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

  !% An implementation of heated dark matter halo profiles.

  !# <darkMatterProfileDMO name="darkMatterProfileDMOHeated">
  !#  <description>Heated dark matter halo profiles.</description>
  !# </darkMatterProfileDMO>

  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough
  use :: Kind_Numbers                , only : kind_int8

  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOHeated
     !% A dark matter halo profile class implementing heated dark matter halos.
     private
     class           (darkMatterProfileDMOClass    ), pointer :: darkMatterProfileDMO_     => null()
     class           (darkMatterProfileHeatingClass), pointer :: darkMatterProfileHeating_ => null()
     integer         (kind=kind_int8               )          :: lastUniqueID
     double precision                                         :: radiusFinalPrevious                , radiusInitialPrevious
     integer                                                  :: nonAnalyticSolver
   contains
     !@ <objectMethods>
     !@   <object>darkMatterProfileDMOHeated</object>
     !@   <objectMethod>
     !@     <method>calculationReset</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(table)\textgreater} node\arginout</arguments>
     !@     <description>Reset memoized calculations.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusInitial</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless *type(treeNode)\textgreater} node\arginout, \doublezero\ radiusFinal\argin</arguments>
     !@     <description>Return the initial radius corresponding to the given final radius in a heated dark matter halo density profile.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                      heatedDestructor
     procedure :: autoHook                          => heatedAutoHook
     procedure :: calculationReset                  => heatedCalculationReset
     procedure :: radiusInitial                     => heatedRadiusInitial
     procedure :: density                           => heatedDensity
     procedure :: densityLogSlope                   => heatedDensityLogSlope
     procedure :: radiusEnclosingDensity            => heatedRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => heatedRadiusEnclosingMass
     procedure :: radialMoment                      => heatedRadialMoment
     procedure :: enclosedMass                      => heatedEnclosedMass
     procedure :: potential                         => heatedPotential
     procedure :: circularVelocity                  => heatedCircularVelocity
     procedure :: circularVelocityMaximum           => heatedCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => heatedRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => heatedRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => heatedRotationNormalization
     procedure :: energy                            => heatedEnergy
     procedure :: energyGrowthRate                  => heatedEnergyGrowthRate
     procedure :: kSpace                            => heatedKSpace
     procedure :: freefallRadius                    => heatedFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => heatedFreefallRadiusIncreaseRate
  end type darkMatterProfileDMOHeated

  interface darkMatterProfileDMOHeated
     !% Constructors for the {\normalfont \ttfamily heated} dark matter halo profile class.
     module procedure heatedConstructorParameters
     module procedure heatedConstructorInternal
  end interface darkMatterProfileDMOHeated

  ! Global variables used in root solving.
  double precision                                      :: heatedRadiusFinal
  type            (treeNode                  ), pointer :: heatedNode
  type            (darkMatterProfileDMOHeated), pointer :: heatedSelf
  !$omp threadprivate(heatedRadiusFinal,heatedNode,heatedSelf)

contains

  function heatedConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily heated} dark matter halo profile class.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (darkMatterProfileDMOHeated    )                :: self
    type   (inputParameters               ), intent(inout) :: parameters
    class  (darkMatterProfileDMOClass     ), pointer       :: darkMatterProfileDMO_
    class  (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    class  (darkMatterProfileHeatingClass ), pointer       :: darkMatterProfileHeating_
    type   (varying_string                )                :: nonAnalyticSolver

    !# <inputParameter>
    !#   <name>nonAnalyticSolver</name>
    !#   <defaultValue>var_str('fallThrough')</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterProfileDMO"     name="darkMatterProfileDMO_"     source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    !# <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    self=darkMatterProfileDMOHeated(enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterProfileDMO_"    />
    !# <objectDestructor name="darkMatterHaloScale_"     />
    !# <objectDestructor name="darkMatterProfileHeating_"/>
    return
  end function heatedConstructorParameters

  function heatedConstructorInternal(nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_) result(self)
    !% Generic constructor for the {\normalfont \ttfamily heated} dark matter profile class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type   (darkMatterProfileDMOHeated   )                        :: self
    class  (darkMatterProfileDMOClass    ), intent(in   ), target :: darkMatterProfileDMO_
    class  (darkMatterHaloScaleClass     ), intent(in   ), target :: darkMatterHaloScale_
    class  (darkMatterProfileHeatingClass), intent(in   ), target :: darkMatterProfileHeating_
    integer                               , intent(in   )         :: nonAnalyticSolver
    !# <constructorAssign variables="nonAnalyticSolver, *darkMatterProfileDMO_, *darkMatterHaloScale_, *darkMatterProfileHeating_"/>

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Galacticus_Error_Report('invalid non-analytic solver type'//{introspection:location})
    ! Construct the object.
    self%lastUniqueID       =-1_kind_int8
    self%radiusFinalPrevious=-huge(0.0d0)
    return
  end function heatedConstructorInternal

  subroutine heatedAutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout) :: self

    call calculationResetEvent%attach(self,heatedCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine heatedAutoHook

  subroutine heatedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily heated} dark matter halo profile class.
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMOHeated), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterProfileDMO_"     />
    !# <objectDestructor name="self%darkMatterHaloScale_"      />
    !# <objectDestructor name="self%darkMatterProfileHeating_" />
    call calculationResetEvent%detach(self,heatedCalculationReset)
    return
  end subroutine heatedDestructor

  subroutine heatedCalculationReset(self,node)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    ! Reset calculations for this profile.
    self%lastUniqueID       =node%uniqueID()
    self%radiusFinalPrevious=-huge(0.0d0)
    return
  end subroutine heatedCalculationReset

  double precision function heatedDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius
    double precision                                            :: radiusInitial , massEnclosed, &
         &                                                         densityInitial, jacobian

    radiusInitial      =self                %radiusInitial(node,radius       )
    massEnclosed       =self%darkMatterProfileDMO_%enclosedMass (node,radiusInitial)
    densityInitial     =self%darkMatterProfileDMO_%density      (node,radiusInitial)
    jacobian           =+1.0d0                                                                                                       &
         &              /(                                                                                                           &
         &                +(                                                                                                         &
         &                  +radius                                                                                                  &
         &                  /radiusInitial                                                                                           &
         &                 )                                                                                                     **2 &
         &                +2.0d0                                                                                                     &
         &                *radius                                                                                                **2 &
         &                /gravitationalConstantGalacticus                                                                           &
         &                /massEnclosed                                                                                              &
         &                *(                                                                                                         &
         &                  +self%darkMatterProfileHeating_%specificEnergyGradient(node,self%darkMatterProfileDMO_,radiusInitial)    &
         &                  -4.0d0                                                                                                   &
         &                  *Pi                                                                                                      &
         &                  *radiusInitial                                                                                       **2 &
         &                  *densityInitial                                                                                          &
         &                  *self%darkMatterProfileHeating_%specificEnergy        (node,self%darkMatterProfileDMO_,radiusInitial)    &
         &                  /massEnclosed                                                                                            &
         &                 )                                                                                                         &
         &               )
    heatedDensity=+densityInitial                                                                                                    &
         &               *(                                                                                                          &
         &                 +radiusInitial                                                                                            &
         &                 /radius                                                                                                   &
         &                )                                                                                                      **2 &
         &               *jacobian
    return
  end function heatedDensity

  double precision function heatedDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedDensityLogSlope=self%darkMatterProfileDMO_%densityLogSlope         (node,radius)
    else
       heatedDensityLogSlope=self                      %densityLogSlopeNumerical(node,radius)
    end if
    return
  end function heatedDensityLogSlope

  double precision function heatedRadiusEnclosingDensity(self,node,density)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout), target :: self
    type            (treeNode                  ), intent(inout), target :: node
    double precision                            , intent(in   )         :: density

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity         (node,density)
    else
       heatedRadiusEnclosingDensity=self                      %radiusEnclosingDensityNumerical(node,density)
    end if
    return
  end function heatedRadiusEnclosingDensity

  double precision function heatedRadiusEnclosingMass(self,node,mass)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    use :: Galactic_Structure_Options  , only : radiusLarge
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout), target :: self
    type            (treeNode                  ), intent(inout), target :: node
    double precision                            , intent(in   )         :: mass
    double precision                                                    :: radiusInitial
    double precision                                                    :: energySpecific

    radiusInitial            =self%darkMatterProfileDMO_    %radiusEnclosingMass(node,mass                                    )
    energySpecific           =self%darkMatterProfileHeating_%specificEnergy     (node,self%darkMatterProfileDMO_,radiusInitial)
    heatedRadiusEnclosingMass=+1.0d0                                                      &
         &                    /                                                           &
         &                    (                                                           &
         &                     +1.0d0/radiusInitial                                       &
         &                     -2.0d0/gravitationalConstantGalacticus/mass*energySpecific &
         &                    )
    ! If the radius found is negative, which means the intial shell has expanded to infinity, return the largest radius.
    if (heatedRadiusEnclosingMass < 0.0d0) heatedRadiusEnclosingMass=radiusLarge
    return
  end function heatedRadiusEnclosingMass

  double precision function heatedRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout)           :: self
    type            (treeNode                  ), intent(inout)           :: node
    double precision                            , intent(in   )           :: moment
    double precision                            , intent(in   ), optional :: radiusMinimum, radiusMaximum

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedRadialMoment=self%darkMatterProfileDMO_%radialMoment         (node,moment,radiusMinimum,radiusMaximum)
    else
       heatedRadialMoment=self                      %radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    end if
    return
  end function heatedRadialMoment

  double precision function heatedEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius

    heatedEnclosedMass=self%darkMatterProfileDMO_%enclosedMass(node,self%radiusInitial(node,radius))
    return
  end function heatedEnclosedMass

  double precision function heatedRadiusInitial(self,node,radiusFinal)
    !% Find the initial radius corresponding to the given {\normalfont \ttfamily radiusFinal} in
    !% the heated dark matter profile.
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout), target  :: self
    type            (treeNode                  ), intent(inout), target  :: node
    double precision                            , intent(in   )          :: radiusFinal
    double precision                            , parameter              :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6
    type            (rootFinder                ), save                   :: finder
    double precision                                                     :: factorExpand
    !$omp threadprivate(finder)

    ! If profile is unheated, the initial radius equals the final radius.
    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_)) then
       heatedRadiusInitial=radiusFinal
       return
    end if
    ! Reset calculations if necessary.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Find the initial radius in the unheated profile.
    if (radiusFinal /= self%radiusFinalPrevious) then
       heatedSelf        => self
       heatedNode        => node
       heatedRadiusFinal =  radiusFinal
       if (.not.finder%isInitialized()) then
          call finder%rootFunction(heatedRadiusInitialRoot                  )
          call finder%tolerance   (toleranceAbsolute      ,toleranceRelative)
       end if
       if (self%radiusFinalPrevious <= -huge(0.0d0) .or. radiusFinal < self%radiusInitialPrevious) then
          ! No previous solution is available, or the requested final radius is smaller than the previous initial radius. In this
          ! case, our guess for the initial radius is the final radius, and we expand the range downward to find a solution.
          call finder%rangeExpand (                                                             &
               &                   rangeExpandUpward            =1.01d0                       , &
               &                   rangeExpandDownward          =0.50d0                       , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
          self%radiusInitialPrevious=finder%find(rootGuess=radiusFinal)
       else
          ! Previous solution exists, and the requested final radius is larger than the previous initial radius. Use the previous
          ! initial radius as a guess for the solution, with range expansion in steps determined by the relative values of the
          ! current and previous final radii. If the current final radius is close to the previous final radius this should give a
          ! guess for the initial radius close to the actual solution.
          if (radiusFinal > self%radiusFinalPrevious) then
             factorExpand=     radiusFinal        /self%radiusFinalPrevious
          else
             factorExpand=self%radiusFinalPrevious/     radiusFinal
          end if
          call finder%rangeExpand (                                                             &
               &                   rangeExpandUpward            =1.0d0*factorExpand           , &
               &                   rangeExpandDownward          =1.0d0/factorExpand           , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
          self%radiusInitialPrevious=finder%find(rootGuess=self%radiusInitialPrevious)
       end if
       self%radiusFinalPrevious=radiusFinal
    end if
    heatedRadiusInitial=self%radiusInitialPrevious
    return
  end function heatedRadiusInitial

  double precision function heatedRadiusInitialRoot(radiusInitial)
    !% Root function used in finding initial radii in heated dark matter halo profiles.
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    double precision, intent(in   ) :: radiusInitial
    double precision                :: massEnclosed

    massEnclosed           =+heatedSelf%darkMatterProfileDMO_    %enclosedMass  (heatedNode,                                 radiusInitial)
    heatedRadiusInitialRoot=+heatedSelf%darkMatterProfileHeating_%specificEnergy(heatedNode,heatedSelf%darkMatterProfileDMO_,radiusInitial) &
         &                  +0.5d0                                                                                                          &
         &                  *gravitationalConstantGalacticus                                                                                &
         &                  *massEnclosed                                                                                                   &
         &                  *(                                                                                                              &
         &                    +1.0d0/heatedRadiusFinal                                                                                      &
         &                    -1.0d0/radiusInitial                                                                                          &
         &                   )
    return
  end function heatedRadiusInitialRoot

  double precision function heatedPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout)           :: self
    type            (treeNode                  ), intent(inout), target   :: node
    double precision                            , intent(in   )           :: radius
    integer                                     , intent(  out), optional :: status

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedPotential=self%darkMatterProfileDMO_%potential         (node,radius,status)
    else
       heatedPotential=self                      %potentialNumerical(node,radius,status)
    end if
    return
  end function heatedPotential

  double precision function heatedCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius

    if (radius > 0.0d0) then
       heatedCircularVelocity=sqrt(                                 &
            &                      +gravitationalConstantGalacticus &
            &                      *self%enclosedMass(node,radius)  &
            &                      /                       radius   &
            &                     )
    else
       heatedCircularVelocity=0.0d0
    end if
    return
  end function heatedCircularVelocity

  double precision function heatedCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedCircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum         (node)
    else
       heatedCircularVelocityMaximum=self                      %circularVelocityMaximumNumerical(node)
    end if
    return
  end function heatedCircularVelocityMaximum

  double precision function heatedRadialVelocityDispersion(self,node,radius)
    !% Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius
    double precision                                            :: radiusInitial           , energySpecific, &
         &                                                         velocityDispersionSquare

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedRadialVelocityDispersion=self%darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
    else
       radiusInitial                 = self%radiusInitial                                     (node                           ,radius       )
       energySpecific                = self%darkMatterProfileHeating_%specificEnergy          (node,self%darkMatterProfileDMO_,radiusInitial)
       velocityDispersionSquare      =+self%darkMatterProfileDMO_    %radialVelocityDispersion(node                           ,radiusInitial)**2 &
            &                         -2.0d0/3.0d0*energySpecific
       heatedRadialVelocityDispersion=sqrt(max(0.0d0,velocityDispersionSquare))
    end if
    return
  end function heatedRadialVelocityDispersion

  double precision function heatedRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), pointer :: node
    double precision                            , intent(in   )          :: specificAngularMomentum

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum         (node,specificAngularMomentum)
    else
       heatedRadiusFromSpecificAngularMomentum=self                      %radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    end if
    return
  end function heatedRadiusFromSpecificAngularMomentum

  double precision function heatedRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedRotationNormalization=self%darkMatterProfileDMO_%rotationNormalization         (node)
    else
       heatedRotationNormalization=self                      %rotationNormalizationNumerical(node)
    end if
    return
  end function heatedRotationNormalization

  double precision function heatedEnergy(self,node)
    !% Return the energy of a heated halo density profile.
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedEnergy=self%darkMatterProfileDMO_%energy         (node)
    else
       heatedEnergy=self                      %energyNumerical(node)
    end if
    return
  end function heatedEnergy

  double precision function heatedEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of a heated halo density profile.
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout)         :: self
    type (treeNode                  ), intent(inout), target :: node

   if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedEnergyGrowthRate=self%darkMatterProfileDMO_%energyGrowthRate         (node)
    else
       heatedEnergyGrowthRate=self                      %energyGrowthRateNumerical(node)
    end if
    return
  end function heatedEnergyGrowthRate

  double precision function heatedKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the heated density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$), using the expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout)         :: self
    type            (treeNode                  ), intent(inout), target :: node
    double precision                            , intent(in   )         :: waveNumber

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedKSpace=self%darkMatterProfileDMO_%kSpace         (node,waveNumber)
    else
       heatedKSpace=self                      %kSpaceNumerical(node,waveNumber)
    end if
    return
  end function heatedKSpace

  double precision function heatedFreefallRadius(self,node,time)
    !% Returns the freefall radius in the heated density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: time

   if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedFreefallRadius=self%darkMatterProfileDMO_%freefallRadius         (node,time)
    else
       heatedFreefallRadius=self                      %freefallRadiusNumerical(node,time)
    end if
    return
  end function heatedFreefallRadius

  double precision function heatedFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the heated density profile at the specified {\normalfont
    !% \ttfamily time} (given in Gyr).
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: time

   if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedFreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate         (node,time)
    else
       heatedFreefallRadiusIncreaseRate=self                      %freefallRadiusIncreaseRateNumerical(node,time)
    end if
    return
  end function heatedFreefallRadiusIncreaseRate
