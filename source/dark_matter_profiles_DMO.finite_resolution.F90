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
  An implementation of dark matter halo profiles with finite resolution (to mimic the effects of resolution in N-body
  simulations for example).
  !!}

  use :: Cosmology_Functions         , only : cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOFiniteResolution">
   <description>
    A dark matter profile DMO class which applies a finite resolution to some other DMO class, typically to mimic the effects
    of finite resolution in an N-body simulation. Specifically, the density profile is given by
    \begin{equation}
     \rho(r) = \rho^\prime(r) \left( 1 + \left[ \frac{\Delta x}{r} \right]^2 \right)^{-1/2},
    \end{equation}
    where $\Delta x$ is the larger of the resolution length, {\normalfont \ttfamily [lengthResolution]}, and the radius in the
    original profile enclosing the mass resolution, {\normalfont \ttfamily [massResolution]}.
  
    Note that this choice was constructed to give a constant density core in an NFW density profile. For a density profile, $\rho^\prime(r)$, which
    rises more steeply than $r^{-1}$ as $r \rightarrow 0$ we will still have a cuspy density profile under this model.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOFiniteResolution
     !!{
     A dark matter halo profile class implementing finiteResolution dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_      => null()
     class           (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_        => null()
     integer                                              :: nonAnalyticSolver
     integer         (kind=kind_int8           )          :: lastUniqueID
     logical                                              :: resolutionIsComoving
     double precision                                     :: lengthResolution                    , massResolution      , &
          &                                                  lengthResolutionPrevious            , enclosedMassPrevious, &
          &                                                  enclosedMassRadiusPrevious
   contains
     !![
     <methods>
       <method description="Return the resolution length in physical units." method="lengthResolutionPhysical"/>
       <method description="Reset memoized calculations."                    method="calculationReset"        />
     </methods>
     !!]
     final     ::                                      finiteResolutionDestructor
     procedure :: autoHook                          => finiteResolutionAutoHook
     procedure :: calculationReset                  => finiteResolutionCalculationReset
     procedure :: density                           => finiteResolutionDensity
     procedure :: densityLogSlope                   => finiteResolutionDensityLogSlope
     procedure :: radiusEnclosingDensity            => finiteResolutionRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => finiteResolutionRadiusEnclosingMass
     procedure :: radialMoment                      => finiteResolutionRadialMoment
     procedure :: enclosedMass                      => finiteResolutionEnclosedMass
     procedure :: potential                         => finiteResolutionPotential
     procedure :: circularVelocity                  => finiteResolutionCircularVelocity
     procedure :: radiusCircularVelocityMaximum     => finiteResolutionRadiusCircularVelocityMaximum
     procedure :: circularVelocityMaximum           => finiteResolutionCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => finiteResolutionRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => finiteResolutionRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => finiteResolutionRotationNormalization
     procedure :: energy                            => finiteResolutionEnergy
     procedure :: kSpace                            => finiteResolutionKSpace
     procedure :: freefallRadius                    => finiteResolutionFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => finiteResolutionFreefallRadiusIncreaseRate
     procedure :: lengthResolutionPhysical          => finiteResolutionLengthResolutionPhysical
  end type darkMatterProfileDMOFiniteResolution

  interface darkMatterProfileDMOFiniteResolution
     !!{
     Constructors for the {\normalfont \ttfamily finiteResolution} dark matter halo profile class.
     !!}
     module procedure finiteResolutionConstructorParameters
     module procedure finiteResolutionConstructorInternal
  end interface darkMatterProfileDMOFiniteResolution

  double precision, parameter :: radiusLengthResolutionRatioMaximum=100.0d0

contains

  function finiteResolutionConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily finiteResolution} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOFiniteResolution)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass           ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), pointer       :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass             ), pointer       :: cosmologyFunctions_
    double precision                                                      :: lengthResolution       , massResolution
    type            (varying_string                      )                :: nonAnalyticSolver
    logical                                                               :: resolutionIsComoving

    !![
    <inputParameter>
      <name>lengthResolution</name>
      <source>parameters</source>
      <description>The resolution length, $\Delta x$.</description>
    </inputParameter>
    <inputParameter>
      <name>massResolution</name>
      <source>parameters</source>
      <description>The resolution mass, $\Delta M$.</description>
    </inputParameter>
    <inputParameter>
      <name>resolutionIsComoving</name>
      <source>parameters</source>
      <description>If true, the resolution length is assumed to be fixed in comoving coordinates, otherwise in physical coordinates.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    !!]
    self=darkMatterProfileDMOFiniteResolution(lengthResolution,massResolution,resolutionIsComoving,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterProfileDMO_,darkMatterHaloScale_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="cosmologyFunctions_"  />
    !!]
    return
  end function finiteResolutionConstructorParameters

  function finiteResolutionConstructorInternal(lengthResolution,massResolution,resolutionIsComoving,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_,cosmologyFunctions_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily finiteResolution} dark matter profile class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (darkMatterProfileDMOFiniteResolution)                        :: self
    class           (darkMatterProfileDMOClass           ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass             ), intent(in   ), target :: cosmologyFunctions_
    double precision                                      , intent(in   )         :: lengthResolution       , massResolution
    integer                                               , intent(in   )         :: nonAnalyticSolver
    logical                                               , intent(in   )         :: resolutionIsComoving
    !![
    <constructorAssign variables="lengthResolution, massResolution, resolutionIsComoving, nonAnalyticSolver, *darkMatterProfileDMO_, *darkMatterHaloScale_, *cosmologyFunctions_"/>
    !!]
    
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    self%lastUniqueID              =-1_kind_int8
    self%genericLastUniqueID       =-1_kind_int8
    self%lengthResolutionPrevious  =-huge(0.0d0)
    self%enclosedMassPrevious      =-huge(0.0d0)
    self%enclosedMassRadiusPrevious=-huge(0.0d0)
    return
  end function finiteResolutionConstructorInternal

  subroutine finiteResolutionAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOFiniteResolution), intent(inout) :: self

    call calculationResetEvent%attach(self,finiteResolutionCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine finiteResolutionAutoHook

  subroutine finiteResolutionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily finiteResolution} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOFiniteResolution), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%cosmologyFunctions_"  />
    !!]
    return
  end subroutine finiteResolutionDestructor

  subroutine finiteResolutionCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    self%lastUniqueID              =node%uniqueID()
    self%genericLastUniqueID       =node%uniqueID()
    self%lengthResolutionPrevious  =-huge(0.0d0)
    self%enclosedMassPrevious      =-huge(0.0d0)
    self%enclosedMassRadiusPrevious=-huge(0.0d0)
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine finiteResolutionCalculationReset

  double precision function finiteResolutionDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    double precision                                                      :: lengthResolution

    lengthResolution       =+self                      %lengthResolutionPhysical(node       )
    finiteResolutionDensity=+self%darkMatterProfileDMO_%density                 (node,radius) &
         &                  /sqrt(                                                            &
         &                        +1.0d0                                                      &
         &                        +(                                                          &
         &                          +lengthResolution                                         &
         &                          /radius                                                   &
         &                         )**2                                                       &
         &                       )
    return
  end function finiteResolutionDensity

  double precision function finiteResolutionDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    double precision                                                      :: lengthResolution

    lengthResolution               =+self                      %lengthResolutionPhysical(node       )
    finiteResolutionDensityLogSlope=+self%darkMatterProfileDMO_%densityLogSlope         (node,radius) &
         &                          +(                                                                &
         &                            +  lengthResolution                                             &
         &                            /  radius                                                       &
         &                           )  **2                                                           &
         &                          /(                                                                &
         &                            +1.0d0                                                          &
         &                            +(                                                              &
         &                              +lengthResolution                                             &
         &                              /radius                                                       &
         &                             )**2                                                           &
         &                           )
    return
  end function finiteResolutionDensityLogSlope

  double precision function finiteResolutionRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: density

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity         (node,density)
    else
       finiteResolutionRadiusEnclosingDensity=self                      %radiusEnclosingDensityNumerical(node,density)
    end if
    return
  end function finiteResolutionRadiusEnclosingDensity

  double precision function finiteResolutionRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: mass

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionRadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass         (node,mass)
    else
       finiteResolutionRadiusEnclosingMass=self                      %radiusEnclosingMassNumerical(node,mass)
    end if
      return
  end function finiteResolutionRadiusEnclosingMass

  double precision function finiteResolutionRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the radial moment of the density profile.
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    double precision                                      , intent(in   )           :: moment
    double precision                                      , intent(in   ), optional :: radiusMinimum, radiusMaximum
  
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionRadialMoment=self%darkMatterProfileDMO_%radialMoment         (node,moment,radiusMinimum,radiusMaximum)
    else
       finiteResolutionRadialMoment=self                      %radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    end if
    return
  end function finiteResolutionRadialMoment

  double precision function finiteResolutionEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius

    if (node%uniqueID() /= self%lastUniqueID              ) call self%calculationReset(node)
    if (     radius     /= self%enclosedMassRadiusPrevious) then
       self%enclosedMassRadiusPrevious=radius
       if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough .or. radius > radiusLengthResolutionRatioMaximum*self%lengthResolutionPhysical(node)) then
          self%enclosedMassPrevious=self%darkMatterProfileDMO_%enclosedMass         (node,radius)
       else
          self%enclosedMassPrevious=self                      %enclosedMassNumerical(node,radius)
       end if
    end if
    finiteResolutionEnclosedMass=self%enclosedMassPrevious
    return
  end function finiteResolutionEnclosedMass

  double precision function finiteResolutionPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout)           :: self
    type            (treeNode                            ), intent(inout), target   :: node
    double precision                                      , intent(in   )           :: radius
    integer                                               , intent(  out), optional :: status

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough .or. radius > radiusLengthResolutionRatioMaximum*self%lengthResolutionPhysical(node)) then
       finiteResolutionPotential=self%darkMatterProfileDMO_%potential         (node,radius,status)
    else
       finiteResolutionPotential=self                      %potentialNumerical(node,radius,status)
    end if
    return
  end function finiteResolutionPotential

  double precision function finiteResolutionCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius

    if (radius > 0.0d0) then
       finiteResolutionCircularVelocity=sqrt(                                 &
            &                                +gravitationalConstantGalacticus &
            &                                *self%enclosedMass(node,radius)  &
            &                                /                       radius   &
            &                               )
    else
       finiteResolutionCircularVelocity=0.0d0
    end if
    return
  end function finiteResolutionCircularVelocity

  double precision function finiteResolutionRadiusCircularVelocityMaximum(self,node)
    !!{
    Returns the radius (in Mpc) at which the maximum circular velocity is acheived in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionRadiusCircularVelocityMaximum=self%darkMatterProfileDMO_%radiusCircularVelocityMaximum         (node)
    else
       finiteResolutionRadiusCircularVelocityMaximum=self                      %radiusCircularVelocityMaximumNumerical(node)
    end if
    return
  end function finiteResolutionRadiusCircularVelocityMaximum

  double precision function finiteResolutionCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionCircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum         (node)
    else
       finiteResolutionCircularVelocityMaximum=self                      %circularVelocityMaximumNumerical(node)
    end if
    return
  end function finiteResolutionCircularVelocityMaximum

  double precision function finiteResolutionRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
 
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough .or. radius > radiusLengthResolutionRatioMaximum*self%lengthResolutionPhysical(node)) then
       finiteResolutionRadialVelocityDispersion=self%darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
    else
       finiteResolutionRadialVelocityDispersion=self                      %radialVelocityDispersionNumerical(node,radius)
    end if
    return
  end function finiteResolutionRadialVelocityDispersion

  double precision function finiteResolutionRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: specificAngularMomentum

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum         (node,specificAngularMomentum)
    else
       finiteResolutionRadiusFromSpecificAngularMomentum=self                      %radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    end if
    return
  end function finiteResolutionRadiusFromSpecificAngularMomentum

  double precision function finiteResolutionRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                                      :: radiusVirial

    radiusVirial                         =+self%darkMatterHaloScale_%radiusVirial(node                                                           )
    finiteResolutionRotationNormalization=+self                     %radialMoment(node,moment=2.0d0,radiusMinimum=0.0d0,radiusMaximum=radiusVirial) &
         &                                /self                     %radialMoment(node,moment=3.0d0,radiusMinimum=0.0d0,radiusMaximum=radiusVirial)
    return
  end function finiteResolutionRotationNormalization

  double precision function finiteResolutionEnergy(self,node)
    !!{
    Return the energy of a finiteResolution halo density profile.
    !!}
    implicit none
    class(darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionEnergy=self%darkMatterProfileDMO_%energy         (node)
    else
       finiteResolutionEnergy=self                      %energyNumerical(node)
    end if
    return
  end function finiteResolutionEnergy

  double precision function finiteResolutionKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the finiteResolution density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$), using the expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout)         :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: waveNumber

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionKSpace=self%darkMatterProfileDMO_%kSpace         (node,waveNumber)
    else
       finiteResolutionKSpace=self                      %kSpaceNumerical(node,waveNumber)
    end if
    return
  end function finiteResolutionKSpace

  double precision function finiteResolutionFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the finiteResolution density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: time

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionFreefallRadius=self%darkMatterProfileDMO_%freefallRadius         (node,time)
    else
       finiteResolutionFreefallRadius=self                      %freefallRadiusNumerical(node,time)
    end if
    return
  end function finiteResolutionFreefallRadius

  double precision function finiteResolutionFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the finiteResolution density profile at the specified {\normalfont
    \ttfamily time} (given in Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOFiniteResolution), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: time

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       finiteResolutionFreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate         (node,time)
    else
       finiteResolutionFreefallRadiusIncreaseRate=self                      %freefallRadiusIncreaseRateNumerical(node,time)
    end if
    return
  end function finiteResolutionFreefallRadiusIncreaseRate

  double precision function finiteResolutionLengthResolutionPhysical(self,node)
    !!{
    Return the resolution length in physical units.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    class(nodeComponentBasic                  ), pointer       :: basic

    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    if (self%lengthResolutionPrevious < 0.0d0) then
       self%lengthResolutionPrevious=self%lengthResolution
       if (self%resolutionIsComoving) then
          basic                            =>                                           node %basic()
          self%lengthResolutionPrevious    =  +self%lengthResolutionPrevious                           &
               &                              *self%cosmologyFunctions_%expansionFactor(basic%time ())
       end if
       if (self%massResolution > 0.0d0)                                                                                        &
            & self%lengthResolutionPrevious=max(                                                                               &
            &                                   self                      %lengthResolutionPrevious                          , &
            &                                   self%darkMatterProfileDMO_%radiusEnclosingMass     (node,self%massResolution)  &
            &                                  )
    end if
    finiteResolutionLengthResolutionPhysical=self%lengthResolutionPrevious
    return
  end function finiteResolutionLengthResolutionPhysical
