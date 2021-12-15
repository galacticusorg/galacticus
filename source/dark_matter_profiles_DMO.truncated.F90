!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  An implementation of truncated dark matter halo profiles.
  !!}

  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOTruncated">
   <description>truncated dark matter halo profiles.</description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOTruncated
     !!{
     A dark matter halo profile class implementing truncated dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_              => null()
     double precision                                     :: radiusFractionalTruncateMinimum                           , radiusFractionalTruncateMaximum
     integer                                              :: nonAnalyticSolver
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8          )           :: lastUniqueID
     ! Stored values of computed quantities.
     double precision                                     :: enclosedMassTruncateMinimumPrevious                       , enclosedMassTruncateMaximumPrevious            , &
          &                                                  enclosingMassRadiusPrevious                               , radialVelocityDispersionTruncateMinimumPrevious, &
          &                                                  radialVelocityDispersionTruncateMinimumUntruncatedPrevious
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
       <method description="Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="truncationFunction" />
     </methods>
     !!]
     final     ::                                      truncatedDestructor
     procedure :: autoHook                          => truncatedAutoHook
     procedure :: calculationReset                  => truncatedCalculationReset
     procedure :: density                           => truncatedDensity
     procedure :: densityLogSlope                   => truncatedDensityLogSlope
     procedure :: radiusEnclosingDensity            => truncatedRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => truncatedRadiusEnclosingMass
     procedure :: radialMoment                      => truncatedRadialMoment
     procedure :: enclosedMass                      => truncatedEnclosedMass
     procedure :: potential                         => truncatedPotential
     procedure :: circularVelocity                  => truncatedCircularVelocity
     procedure :: circularVelocityMaximum           => truncatedCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => truncatedRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => truncatedRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => truncatedRotationNormalization
     procedure :: energy                            => truncatedEnergy
     procedure :: kSpace                            => truncatedKSpace
     procedure :: freefallRadius                    => truncatedFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => truncatedFreefallRadiusIncreaseRate
     procedure :: truncationFunction                => truncatedTruncationFunction
  end type darkMatterProfileDMOTruncated

  interface darkMatterProfileDMOTruncated
     !!{
     Constructors for the {\normalfont \ttfamily truncated} dark matter halo profile class.
     !!}
     module procedure truncatedConstructorParameters
     module procedure truncatedConstructorInternal
  end interface darkMatterProfileDMOTruncated

contains

  function truncatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily truncated} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOTruncated)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass    ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
    type            (varying_string               )                :: nonAnalyticSolver
    double precision                                               :: radiusFractionalTruncateMinimum, radiusFractionalTruncateMaximum

    !![
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusFractionalTruncateMinimum</name>
      <defaultValue>2.0d0</defaultValue>
      <source>parameters</source>
      <description>The minimum radius (in units of the virial radius) to begin truncating the density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusFractionalTruncateMaximum</name>
      <defaultValue>4.0d0</defaultValue>
      <source>parameters</source>
      <description>The maximum radius (in units of the virial radius) to finish truncating the density profile.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO"   name="darkMatterProfileDMO_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOTruncated(radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"  />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function truncatedConstructorParameters

  function truncatedConstructorInternal(radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily truncated} dark matter profile class.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (darkMatterProfileDMOTruncated)                        :: self
    class           (darkMatterProfileDMOClass    ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass     ), intent(in   ), target :: darkMatterHaloScale_
    double precision                               , intent(in   )         :: radiusFractionalTruncateMinimum, radiusFractionalTruncateMaximum
    integer                                        , intent(in   )         :: nonAnalyticSolver
    !![
    <constructorAssign variables="radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,nonAnalyticSolver,*darkMatterProfileDMO_,*darkMatterHaloScale_"/>
    !!]

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Galacticus_Error_Report('invalid non-analytic solver type'//{introspection:location})
    self%lastUniqueID       =-1_kind_int8
    self%genericLastUniqueID=-1_kind_int8
    return
  end function truncatedConstructorInternal

  subroutine truncatedAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout) :: self

    call calculationResetEvent%attach(self,truncatedCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine truncatedAutoHook

  subroutine truncatedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily truncated} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMOTruncated), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    if (calculationResetEvent%isAttached(self,truncatedCalculationReset)) call calculationResetEvent%detach(self,truncatedCalculationReset)
    return
  end subroutine truncatedDestructor

  subroutine truncatedCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node

    self%lastUniqueID                                              =node%uniqueID()
    self%genericLastUniqueID                                       =node%uniqueID()
    self%enclosingMassRadiusPrevious                               =-1.0d0
    self%enclosedMassTruncateMinimumPrevious                       =-1.0d0
    self%enclosedMassTruncateMaximumPrevious                       =-1.0d0
    self%radialVelocityDispersionTruncateMinimumPrevious           =-1.0d0
    self%radialVelocityDispersionTruncateMinimumUntruncatedPrevious=-1.0d0
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine truncatedCalculationReset

  subroutine truncatedTruncationFunction(self,node,radius,x,multiplier,multiplierGradient)
    !!{
    Return the scaled truncation radial coordinate, and the truncation multiplier.
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout)           :: self
    type            (treeNode                     ), intent(inout)           :: node
    double precision                               , intent(in   )           :: radius
    double precision                               , intent(  out), optional :: x                 , multiplier, &
         &                                                                      multiplierGradient
    double precision                                                         :: radiusVirial      , x_

    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    radiusVirial=self%darkMatterHaloScale_%radiusVirial(node)
    if      (radius <= radiusVirial*self%radiusFractionalTruncateMinimum) then
       if (present(x                 )) x                 =+0.0d0
       if (present(multiplier        )) multiplier        =+1.0d0
       if (present(multiplierGradient)) multiplierGradient=+0.0d0
    else if (radius >= radiusVirial*self%radiusFractionalTruncateMaximum) then
       if (present(x                 )) x                 =+1.0d0
       if (present(multiplier        )) multiplier        =+0.0d0
       if (present(multiplierGradient)) multiplierGradient=+0.0d0
    else
       x_                                                 =+(     radius                         /radiusVirial-self%radiusFractionalTruncateMinimum) &
            &                                              /(self%radiusFractionalTruncateMaximum             -self%radiusFractionalTruncateMinimum)
       if (present(x                 )) x                 =x_
       if (present(multiplier        )) multiplier        =  +1.0d0       &
            &                                                -3.0d0*x_**2 &
            &                                                +2.0d0*x_**3
       if (present(multiplierGradient)) multiplierGradient=+(                                                                                        &
            &                                                -6.0d0*x_                                                                               &
            &                                                +6.0d0*x_**2                                                                            &
            &                                               )                                                                                        &
            &                                              /                                      radiusvirial                                       &
            &                                              /(self%radiusFractionalTruncateMaximum             -self%radiusFractionalTruncateMinimum)
    end if
    return
  end subroutine truncatedTruncationFunction

  double precision function truncatedDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    double precision                                               :: multiplier

    call self%truncationFunction(node,radius,multiplier=multiplier)
    truncatedDensity=+self%darkMatterProfileDMO_%density(node,radius) &
         &           *multiplier
    return
  end function truncatedDensity

  double precision function truncatedDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    double precision                                               :: multiplier, multiplierGradient

    call self%truncationFunction(node,radius,multiplier=multiplier,multiplierGradient=multiplierGradient)
    if (multiplier > 0.0d0) then
       truncatedDensityLogSlope=+self%darkMatterProfileDMO_%densityLogSlope(node,radius) &
            &                   +radius                                                  &
            &                   *multiplierGradient                                      &
            &                   /multiplier
    else
       truncatedDensityLogSlope=+0.0d0
    end if
    return
  end function truncatedDensityLogSlope

  double precision function truncatedRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout), target :: self
    type            (treeNode                     ), intent(inout), target :: node
    double precision                               , intent(in   )         :: density

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity         (node,density)
    else
       truncatedRadiusEnclosingDensity=self                      %radiusEnclosingDensityNumerical(node,density)
    end if
    return
  end function truncatedRadiusEnclosingDensity

  double precision function truncatedRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout), target :: self
    type            (treeNode                     ), intent(inout), target :: node
    double precision                               , intent(in   )         :: mass
    double precision                                                       :: radiusVirial, radiusTruncateMinimum

    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    radiusVirial         =self%darkMatterHaloScale_%radiusVirial(node)
    radiusTruncateMinimum=radiusVirial*self%radiusFractionalTruncateMinimum
    if (self%enclosedMassTruncateMinimumPrevious < 0.0d0) then
       self%enclosedMassTruncateMinimumPrevious=self%enclosedMass(node,radiusTruncateMinimum)
    end if
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough .or. mass <= self%enclosedMassTruncateMinimumPrevious) then
       truncatedRadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass         (node,mass)
    else
       truncatedRadiusEnclosingMass=self                      %radiusEnclosingMassNumerical(node,mass)
    end if
    return
  end function truncatedRadiusEnclosingMass

  double precision function truncatedRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout)           :: self
    type            (treeNode                     ), intent(inout)           :: node
    double precision                               , intent(in   )           :: moment
    double precision                               , intent(in   ), optional :: radiusMinimum, radiusMaximum

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedRadialMoment=self%darkMatterProfileDMO_%radialMoment         (node,moment,radiusMinimum,radiusMaximum)
    else
       truncatedRadialMoment=self                      %radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    end if
    return
  end function truncatedRadialMoment

  double precision function truncatedEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    double precision                                               :: radiusVirial, radiusTruncateMinimum

    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    radiusVirial         =self%darkMatterHaloScale_%radiusVirial(node)
    radiusTruncateMinimum=radiusVirial*self%radiusFractionalTruncateMinimum
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough .or. radius <= radiusTruncateMinimum) then
       truncatedEnclosedMass=+self%darkMatterProfileDMO_%enclosedMass                   (node,radius                      )
    else
       truncatedEnclosedMass=+self%darkMatterProfileDMO_%enclosedMass                   (node,radiusTruncateMinimum       ) &
            &                +self                      %enclosedMassDifferenceNumerical(node,radiusTruncateMinimum,radius)
    end if
    return
  end function truncatedEnclosedMass

  double precision function truncatedPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout)           :: self
    type            (treeNode                     ), intent(inout), target   :: node
    double precision                               , intent(in   )           :: radius
    integer                                        , intent(  out), optional :: status

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedPotential=self%darkMatterProfileDMO_%potential         (node,radius,status)
    else
       truncatedPotential=self                      %potentialNumerical(node,radius,status)
    end if
    return
  end function truncatedPotential

  double precision function truncatedCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedCircularVelocity=self%darkMatterProfileDMO_%circularVelocity         (node,radius)
    else
       truncatedCircularVelocity=self                      %circularVelocityNumerical(node,radius)
    end if
    return
  end function truncatedCircularVelocity

  double precision function truncatedCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedCircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum         (node)
    else
       truncatedCircularVelocityMaximum=self                      %circularVelocityMaximumNumerical(node)
    end if
    return
  end function truncatedCircularVelocityMaximum

  double precision function truncatedRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    double precision                                               :: radiusVirial, radiusTruncateMinimum

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedRadialVelocityDispersion=self%darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
    else
       if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
       radiusVirial         =self%darkMatterHaloScale_%radiusVirial(node)
       radiusTruncateMinimum=radiusVirial*self%radiusFractionalTruncateMinimum
       if (radius >= radiusTruncateMinimum) then
          truncatedRadialVelocityDispersion=self%radialVelocityDispersionNumerical(node,radius)
       else
          if (self%radialVelocityDispersionTruncateMinimumPrevious < 0.0d0 .or. self%radialVelocityDispersionTruncateMinimumUntruncatedPrevious < 0.0d0) then
             self%radialVelocityDispersionTruncateMinimumPrevious           =self%radialVelocityDispersionNumerical             (node,radiusTruncateMinimum)
             self%radialVelocityDispersionTruncateMinimumUntruncatedPrevious=self%darkMatterProfileDMO_%radialVelocityDispersion(node,radiusTruncateMinimum)
          end if
          truncatedRadialVelocityDispersion=sqrt(                                                                                       &
               &                                 +self%darkMatterProfileDMO_%radialVelocityDispersion   (node,radius               )**2 &
               &                                 +self%density                                          (node,radiusTruncateMinimum)    &
               &                                 /self%density                                          (node,radius               )    &
               &                                 *(                                                                                     &
               &                                   +self%radialVelocityDispersionTruncateMinimumPrevious                            **2 &
               &                                   -self%radialVelocityDispersionTruncateMinimumUntruncatedPrevious                 **2 &
               &                                  )                                                                                     &
               &                                )
       end if
    end if
    return
  end function truncatedRadialVelocityDispersion

  double precision function truncatedRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: specificAngularMomentum

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum         (node,specificAngularMomentum)
    else
       truncatedRadiusFromSpecificAngularMomentum=self                      %radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    end if
    return
  end function truncatedRadiusFromSpecificAngularMomentum

  double precision function truncatedRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedRotationNormalization=self%darkMatterProfileDMO_%rotationNormalization         (node)
    else
       truncatedRotationNormalization=self                      %rotationNormalizationNumerical(node)
    end if
    return
  end function truncatedRotationNormalization

  double precision function truncatedEnergy(self,node)
    !!{
    Return the energy of a truncated halo density profile.
    !!}
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedEnergy=self%darkMatterProfileDMO_%energy         (node)
    else
       truncatedEnergy=self                      %energyNumerical(node)
    end if
    return
  end function truncatedEnergy

  double precision function truncatedKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the truncated density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout)         :: self
    type            (treeNode                     ), intent(inout), target :: node
    double precision                               , intent(in   )         :: waveNumber

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedKSpace=self%darkMatterProfileDMO_%kSpace         (node,waveNumber)
    else
       truncatedKSpace=self                      %kSpaceNumerical(node,waveNumber)
    end if
    return
  end function truncatedKSpace

  double precision function truncatedFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the truncated density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout), target :: self
    type            (treeNode                     ), intent(inout), target :: node
    double precision                               , intent(in   )         :: time

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedFreefallRadius=self%darkMatterProfileDMO_%freefallRadius         (node,time)
    else
       truncatedFreefallRadius=self                      %freefallRadiusNumerical(node,time)
    end if
    return
  end function truncatedFreefallRadius

  double precision function truncatedFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the truncated density profile at the specified {\normalfont
    \ttfamily time} (given in Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout), target :: self
    type            (treeNode                     ), intent(inout), target :: node
    double precision                               , intent(in   )         :: time

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedFreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate         (node,time)
    else
       truncatedFreefallRadiusIncreaseRate=self                      %freefallRadiusIncreaseRateNumerical(node,time)
    end if
    return
  end function truncatedFreefallRadiusIncreaseRate
