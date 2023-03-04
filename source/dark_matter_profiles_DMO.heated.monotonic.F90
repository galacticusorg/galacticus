!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  An implementation of heated dark matter halo profiles based on the energy ordering of shells.
  !!}

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOHeatedMonotonic">
   <description>
     A dark matter profile DMO class in which dark matter halos start out with a density profile defined by another {\normalfont
     \ttfamily darkMatterProfileDMO}. This profile is then modified by heating, under the assumption that the
     energy of a shell of mass before and after heating are related by
     \begin{equation}
     -{ \mathrm{G} M^\prime(r^\prime) \over r^\prime } = -{ \mathrm{G} M(r) \over r } + 2 \epsilon(r),
     \end{equation}    
     where $M(r)$ is the mass enclosed within a radius $r$, and $\epsilon(r)$ represents the specific heating in the shell
     initially at radius $r$. Primes indicate values after heating, while unprimed variables indicate quantities prior to
     heating.

     The above equation can be re-written as
     \begin{equation}
     -r^{\prime -1} = -r^{-1} + \xi(r),
     \end{equation}     
     where $\xi(r) = 2 \epsilon(r)/[\mathrm{G} M(r)/r]$ measures the perturbation to the shell. To avoid shell crossing a
     monotonicity relation $r_1 &lt; r_2 \implies \xi(r_1) \le \xi(r_2)$ is enforced by starting at large radius and stepping inward,
     enforcing the condition in the next innermost shell as necessary.
     
     Not all methods have analytic solutions for this profile. If {\normalfont \ttfamily [nonAnalyticSolver]}$=${\normalfont
     \ttfamily fallThrough} then attempts to call these methods in heated profiles will simply return the result from the
     unheated profile, otherwise a numerical calculation is performed.
   </description>
  </darkMatterProfileDMO>
  !!]

  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversType, enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough
  use :: Kind_Numbers                , only : kind_int8
  use :: Numerical_Interpolation     , only : interpolator

  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOHeatedMonotonic
     !!{
     A dark matter halo profile class implementing heated dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass        ), pointer     :: darkMatterProfileDMO_     => null()
     class           (darkMatterProfileHeatingClass    ), pointer     :: darkMatterProfileHeating_ => null()
     integer         (kind=kind_int8                   )              :: lastUniqueID
     type            (enumerationNonAnalyticSolversType)              :: nonAnalyticSolver
     double precision                                                 :: radiusInitialMinimum               , radiusInitialMaximum, &
          &                                                              radiusFinalMinimum                 , radiusFinalMaximum
     type            (interpolator                     ), allocatable :: massProfile
     logical                                                          :: isBound
   contains
     !![
     <methods>
       <method description="Reset memoized calculations."               method="calculationReset"/>
       <method description="Compute a solution for the heated profile." method="computeSolution" />
     </methods>
     !!]
     final     ::                                      heatedMonotonicDestructor
     procedure :: autoHook                          => heatedMonotonicAutoHook
     procedure :: calculationReset                  => heatedMonotonicCalculationReset
     procedure :: computeSolution                   => heatedMonotonicComputeSolution
     procedure :: density                           => heatedMonotonicDensity
     procedure :: densityLogSlope                   => heatedMonotonicDensityLogSlope
     procedure :: radiusEnclosingDensity            => heatedMonotonicRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => heatedMonotonicRadiusEnclosingMass
     procedure :: radialMoment                      => heatedMonotonicRadialMoment
     procedure :: enclosedMass                      => heatedMonotonicEnclosedMass
     procedure :: potential                         => heatedMonotonicPotential
     procedure :: circularVelocity                  => heatedMonotonicCircularVelocity
     procedure :: radiusCircularVelocityMaximum     => heatedMonotonicRadiusCircularVelocityMaximum
     procedure :: circularVelocityMaximum           => heatedMonotonicCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => heatedMonotonicRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => heatedMonotonicRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => heatedMonotonicRotationNormalization
     procedure :: energy                            => heatedMonotonicEnergy
     procedure :: kSpace                            => heatedMonotonicKSpace
     procedure :: freefallRadius                    => heatedMonotonicFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => heatedMonotonicFreefallRadiusIncreaseRate
  end type darkMatterProfileDMOHeatedMonotonic

  interface darkMatterProfileDMOHeatedMonotonic
     !!{
     Constructors for the {\normalfont \ttfamily heated} dark matter halo profile class.
     !!}
     module procedure heatedMonotonicConstructorParameters
     module procedure heatedMonotonicConstructorInternal
  end interface darkMatterProfileDMOHeatedMonotonic

contains

  function heatedMonotonicConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily heatedMonotonic} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOHeatedMonotonic)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass          ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileHeatingClass      ), pointer       :: darkMatterProfileHeating_
    type            (varying_string                     )                :: nonAnalyticSolver
    double precision                                                     :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, &
         &                                                                  toleranceRelativePotential

    !![
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersion</name>
      <defaultValue>1.0d-6</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersionMaximum</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum allowed relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles before aborting.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativePotential</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum allowed relative tolerance to use in numerical solutions for the gravitational potential in dark-matter-only density profiles before aborting.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO"     name="darkMatterProfileDMO_"     source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOHeatedMonotonic(enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,toleranceRelativePotential,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"    />
    <objectDestructor name="darkMatterHaloScale_"     />
    <objectDestructor name="darkMatterProfileHeating_"/>
    !!]
    return
  end function heatedMonotonicConstructorParameters

  function heatedMonotonicConstructorInternal(nonAnalyticSolver,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,toleranceRelativePotential,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily heatedMonotonic} dark matter profile class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (darkMatterProfileDMOHeatedMonotonic)                        :: self
    class           (darkMatterProfileDMOClass          ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileHeatingClass      ), intent(in   ), target :: darkMatterProfileHeating_
    type            (enumerationNonAnalyticSolversType  ), intent(in   )         :: nonAnalyticSolver
    double precision                                     , intent(in   )         :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, &
         &                                                                          toleranceRelativePotential
    !![
    <constructorAssign variables="nonAnalyticSolver, toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, toleranceRelativePotential, *darkMatterProfileDMO_, *darkMatterHaloScale_, *darkMatterProfileHeating_"/>
    !!]

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    ! Construct the object.
    self%lastUniqueID        =-1_kind_int8
    self%genericLastUniqueID =-1_kind_int8
    self%radiusInitialMinimum=+huge(0.0d0)
    self%radiusInitialMaximum=-huge(0.0d0)
    self%radiusFinalMinimum  =+huge(0.0d0)
    self%radiusFinalMaximum  =-huge(0.0d0)
    self%isBound             =.true.
    return
  end function heatedMonotonicConstructorInternal

  subroutine heatedMonotonicAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self

    call calculationResetEvent%attach(self,heatedMonotonicCalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileDMOHeatedMonotonic')
    return
  end subroutine heatedMonotonicAutoHook

  subroutine heatedMonotonicDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily heatedMonotonic} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"     />
    <objectDestructor name="self%darkMatterHaloScale_"      />
    <objectDestructor name="self%darkMatterProfileHeating_" />
    !!]
    if (calculationResetEvent%isAttached(self,heatedMonotonicCalculationReset)) call calculationResetEvent%detach(self,heatedMonotonicCalculationReset)
    return
  end subroutine heatedMonotonicDestructor

  subroutine heatedMonotonicCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node

    ! Reset calculations for this profile.
    self%lastUniqueID                                =node%uniqueID()
    self%genericLastUniqueID                         =node%uniqueID() 
    self%isBound                                     =.true.
    self%radiusInitialMinimum                        =+huge(0.0d0)
    self%radiusInitialMaximum                        =-huge(0.0d0)
    self%radiusFinalMinimum                          =+huge(0.0d0)
    self%radiusFinalMaximum                          =-huge(0.0d0)
    self%genericEnclosedMassRadiusMinimum            =+huge(0.0d0)
    self%genericEnclosedMassRadiusMaximum            =-huge(0.0d0)
    self%genericEnclosedMassRadiusMinimum            =+huge(0.0d0)
    self%genericEnclosedMassRadiusMaximum            =-huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMinimum=+huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMaximum=-huge(0.0d0)
    if (allocated(self%massProfile                            )) deallocate(self%massProfile                            )
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine heatedMonotonicCalculationReset

  subroutine heatedMonotonicComputeSolution(self,node,radius)
    !!{
    Compute the solution for the heated density profile.
    !!}
    use :: Numerical_Ranges                , only : Make_Range                     , rangeTypeLogarithmic
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Table_Labels                    , only : extrapolationTypeFix           , extrapolationTypeZero
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout)             :: self
    type            (treeNode                           ), intent(inout)             :: node
    double precision                                     , intent(in   )             :: radius
    double precision                                     , parameter                 :: radiusFractionMinimum=1.0d-6, radiusFractionMaximum=10.0d0
    integer                                              , parameter                 :: countPerDecadeRadius =100
    double precision                                     , allocatable, dimension(:) :: massEnclosed                , massShell                   , &
         &                                                                              radiusInitial               , radiusFinal                 , &
         &                                                                              energyFinal                 , perturbation
    logical                                              , allocatable, dimension(:) :: isBound
    integer                                                                          :: i                           , countRadii

    ! Determine if we need to retabulate.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Nothing to do if profile is already tabulated.
    if (allocated(self%massProfile)) return
    ! Choose extent of radii at which to tabulate the initial profile.
    self%radiusInitialMinimum=radiusFractionMinimum*self%darkMatterHaloScale_%radiusVirial(node)
    self%radiusInitialMaximum=radiusFractionMaximum*self%darkMatterHaloScale_%radiusVirial(node)
    ! Build grid of radii.
    countRadii=int(log10(self%radiusInitialMaximum/self%radiusInitialMinimum)*dble(countPerDecadeRadius)+1.0d0)
    if (allocated(radiusInitial)) then
       deallocate(radiusInitial)
       deallocate(radiusFinal  )
       deallocate(massEnclosed )
       deallocate(massShell    )
       deallocate(energyFinal  )
       deallocate(perturbation )
    end if
    allocate(radiusInitial(countRadii))
    allocate(radiusFinal  (countRadii))
    allocate(massEnclosed (countRadii))
    allocate(massShell    (countRadii))
    allocate(energyFinal  (countRadii))
    allocate(perturbation (countRadii))
    radiusInitial=Make_Range(self%radiusInitialMinimum,self%radiusInitialMaximum,countRadii,rangeTypeLogarithmic)
    ! Evaluate masses and energies of shells.
    do i=countRadii,1,-1
       massEnclosed(i)=+self%darkMatterProfileDMO_     %enclosedMass  (node,radiusInitial(i)                           )
       perturbation(i)=+2.0d0                                                                                            &
            &          *self%darkMatterProfileHeating_ %specificEnergy(node,radiusInitial(i),self%darkMatterProfileDMO_) &
            &          /gravitationalConstantGalacticus                                                                  &
            &          /                                                    massEnclosed (i)                             &
            &          *                                                    radiusInitial(i)
       ! Limit the perturbation to avoid shell-crossing.
       if (i < countRadii)                              &
            & perturbation(i)=min(                      &
            &                     +perturbation  (i  ), &
            &                     +1.0d0                &
            &                     -radiusInitial (i  )  &
            &                     /radiusInitial (i+1)  &
            &                     *(                    &
            &                       +massEnclosed(i  )  &
            &                       /massEnclosed(i+1)  &
            &                      )**(-1.0d0/3.0d0)    &
            &                     *(                    &
            &                       +1.0d0              &
            &                       -perturbation(i+1)  &
            &                      )                    &
            &                    )
    end do
    ! Compute the final energy of the heated profile.
    energyFinal=+gravitationalConstantGalacticus &
         &      *massEnclosed                    &
         &      /radiusInitial                   &
         &      *(                               &
         &        -1.0d0                         &
         &        +perturbation                  &
         &       )    
    ! Find shell masses.
    massShell(1           )=+massEnclosed(1             )
    massShell(2:countRadii)=+massEnclosed(2:countRadii  ) &
         &                  -massEnclosed(1:countRadii-1)
    ! Evaluation boundedness.
    isBound= energyFinal < 0.0d0 &
         &  .and.                &
         &   massShell   > 0.0d0
    ! Find final radii.
    where (isBound)
       radiusFinal=-gravitationalConstantGalacticus &
            &      *massEnclosed                    &
            &      /energyFinal
    elsewhere
       radiusFinal=+huge(0.0d0)
    end where
    ! Build the final profile interpolator.
    self%isBound=count(isBound) > 2
    if (self%isBound) then
       self%radiusFinalMinimum =minval(radiusFinal ,mask=isBound)
       self%radiusFinalMaximum =maxval(radiusFinal ,mask=isBound)
       ! Construct the interpolator.
       if (allocated(self%massProfile)) deallocate(self%massProfile)
       allocate(self%massProfile)
       self%massProfile=interpolator(                                                        &
            &                        x                =log(pack(radiusFinal ,mask=isBound)), &
            &                        y                =log(pack(massEnclosed,mask=isBound)), &
            &                        extrapolationType=extrapolationTypeFix                  &
            &                       )
    end if
    return
  end subroutine heatedMonotonicComputeSolution

  double precision function heatedMonotonicEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: radius

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_)) then
       ! No heating - use the unheated solution.
       heatedMonotonicEnclosedMass=self%darkMatterProfileDMO_%enclosedMass(node,radius)
    else if (radius <= 0.0d0) then
       ! Non-positive radius, mass must be zero.
       heatedMonotonicEnclosedMass=0.0d0
    else
       ! Compute the solution (as needed).
       call self%computeSolution(node,radius)
       ! For bound halos, interpolate to find the enclosed mass. For unbound halos the enclosed mass is zero.
       if (self%isBound) then
          if (radius < self%radiusFinalMinimum) then
             ! Assume constant density below the minimum radius.
             heatedMonotonicEnclosedMass=+exp(self%massProfile%interpolate(log(self%radiusFinalMinimum))) &
                  &                      *(                                                               &
                  &                        +     radius                                                   &
                  &                        /self%radiusFinalMinimum                                       &
                  &                       )**3
          else
             heatedMonotonicEnclosedMass=+exp(self%massProfile%interpolate(log(radius)))
          end if
       else
          heatedMonotonicEnclosedMass=0.0d0
       end if
    end if
    return
  end function heatedMonotonicEnclosedMass

  double precision function heatedMonotonicDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: radius
    double precision                                                     :: radius_
    
    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_)) then
       ! No heating - use the unheated solution.
       heatedMonotonicDensity=self%darkMatterProfileDMO_%density(node,radius)
    else
       call self%computeSolution(node,radius)
       ! For bound halos, interpolate to find the density. For unbound halos the density is zero.
       if (self%isBound) then
          radius_               =max(radius,self%radiusFinalMinimum)
          heatedMonotonicDensity=+   self%massProfile%derivative (log(radius_))  &
               &                *exp(self%massProfile%interpolate(log(radius_))) &
               &                /4.0d0                                           &
               &                /Pi                                              &
               &                /radius**3
       else
          heatedMonotonicDensity=+0.0d0
       end if
    end if
    return
  end function heatedMonotonicDensity

  double precision function heatedMonotonicDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: radius

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicDensityLogSlope=self%darkMatterProfileDMO_%densityLogSlope         (node,radius)
    else
       heatedMonotonicDensityLogSlope=self                      %densityLogSlopeNumerical(node,radius)
    end if
    return
  end function heatedMonotonicDensityLogSlope

  double precision function heatedMonotonicRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout), target :: self
    type            (treeNode                           ), intent(inout), target :: node
    double precision                                     , intent(in   )         :: density

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity         (node,density)
    else
       heatedMonotonicRadiusEnclosingDensity=self                      %radiusEnclosingDensityNumerical(node,density)
    end if
    return
  end function heatedMonotonicRadiusEnclosingDensity

  double precision function heatedMonotonicRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout), target :: self
    type            (treeNode                           ), intent(inout), target :: node
    double precision                                     , intent(in   )         :: mass

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicRadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass         (node,mass)
    else
       heatedMonotonicRadiusEnclosingMass=self                      %radiusEnclosingMassNumerical(node,mass)
    end if
    return
  end function heatedMonotonicRadiusEnclosingMass

  double precision function heatedMonotonicRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout)           :: self
    type            (treeNode                           ), intent(inout)           :: node
    double precision                                     , intent(in   )           :: moment
    double precision                                     , intent(in   ), optional :: radiusMinimum, radiusMaximum

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicRadialMoment=self%darkMatterProfileDMO_%radialMoment         (node,moment,radiusMinimum,radiusMaximum)
    else
       heatedMonotonicRadialMoment=self                      %radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    end if
    return
  end function heatedMonotonicRadialMoment

  double precision function heatedMonotonicPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout)           :: self
    type            (treeNode                           ), intent(inout), target   :: node
    double precision                                     , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType  ), intent(  out), optional :: status

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicPotential=self%darkMatterProfileDMO_%potential         (node,radius,status)
    else
       heatedMonotonicPotential=self                      %potentialNumerical(node,radius,status)
    end if
    return
  end function heatedMonotonicPotential

  double precision function heatedMonotonicCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: radius

    if (radius > 0.0d0) then
       heatedMonotonicCircularVelocity=sqrt(                                 &
            &                              +gravitationalConstantGalacticus &
            &                              *self%enclosedMass(node,radius)  &
            &                              /                       radius   &
            &                             )
    else
       heatedMonotonicCircularVelocity=0.0d0
    end if
    return
  end function heatedMonotonicCircularVelocity

  double precision function heatedMonotonicRadiusCircularVelocityMaximum(self,node)
    !!{
    Returns the radius (in Mpc) at which the maximum circular velocity is achieved in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicRadiusCircularVelocityMaximum=self%darkMatterProfileDMO_%radiusCircularVelocityMaximum         (node)
    else
       heatedMonotonicRadiusCircularVelocityMaximum=self                      %radiusCircularVelocityMaximumNumerical(node)
    end if
    return
  end function heatedMonotonicRadiusCircularVelocityMaximum

  double precision function heatedMonotonicCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicCircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum         (node)
    else
       heatedMonotonicCircularVelocityMaximum=self                      %circularVelocityMaximumNumerical(node)
    end if
    return
  end function heatedMonotonicCircularVelocityMaximum

  double precision function heatedMonotonicRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: radius

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicRadialVelocityDispersion=self%darkMatterProfileDMO_%radialVelocityDispersion         (node,radius)
    else
       heatedMonotonicRadialVelocityDispersion=self                      %radialVelocityDispersionNumerical(node,radius)
    end if
    return
  end function heatedMonotonicRadialVelocityDispersion

  double precision function heatedMonotonicRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: specificAngularMomentum

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum         (node,specificAngularMomentum)
    else
       heatedMonotonicRadiusFromSpecificAngularMomentum=self                      %radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    end if
    return
  end function heatedMonotonicRadiusFromSpecificAngularMomentum

  double precision function heatedMonotonicRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class(darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicRotationNormalization=self%darkMatterProfileDMO_%rotationNormalization         (node)
    else
       heatedMonotonicRotationNormalization=self                      %rotationNormalizationNumerical(node)
    end if
    return
  end function heatedMonotonicRotationNormalization

  double precision function heatedMonotonicEnergy(self,node)
    !!{
    Return the energy of a heated halo density profile.
    !!}
    implicit none
    class(darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicEnergy=self%darkMatterProfileDMO_%energy         (node)
    else
       heatedMonotonicEnergy=self                      %energyNumerical(node)
    end if
    return
  end function heatedMonotonicEnergy

  double precision function heatedMonotonicKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the heated density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$), using the expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout)         :: self
    type            (treeNode                           ), intent(inout), target :: node
    double precision                                     , intent(in   )         :: waveNumber

    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicKSpace=self%darkMatterProfileDMO_%kSpace         (node,waveNumber)
    else
       heatedMonotonicKSpace=self                      %kSpaceNumerical(node,waveNumber)
    end if
    return
  end function heatedMonotonicKSpace

  double precision function heatedMonotonicFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the heated density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout), target :: self
    type            (treeNode                           ), intent(inout), target :: node
    double precision                                     , intent(in   )         :: time

   if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicFreefallRadius=self%darkMatterProfileDMO_%freefallRadius         (node,time)
    else
       heatedMonotonicFreefallRadius=self                      %freefallRadiusNumerical(node,time)
    end if
    return
  end function heatedMonotonicFreefallRadius

  double precision function heatedMonotonicFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the heated density profile at the specified {\normalfont
    \ttfamily time} (given in Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout), target :: self
    type            (treeNode                           ), intent(inout), target :: node
    double precision                                     , intent(in   )         :: time

   if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_) .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       heatedMonotonicFreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate         (node,time)
    else
       heatedMonotonicFreefallRadiusIncreaseRate=self                      %freefallRadiusIncreaseRateNumerical(node,time)
    end if
    return
  end function heatedMonotonicFreefallRadiusIncreaseRate
