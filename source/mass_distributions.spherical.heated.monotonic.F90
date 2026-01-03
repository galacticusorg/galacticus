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
  An implementation of heated dark matter halo profiles based on the energy ordering of shells.
  !!}

  use :: Numerical_Interpolation, only : interpolator

  !![
  <massDistribution name="massDistributionSphericalHeatedMonotonic">
   <description>
     A mass distribution class in which dark matter halos start out with a density profile defined by another {\normalfont
     \ttfamily massDistributionClass}. This profile is then modified by heating, under the assumption that the
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
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalDecorator) :: massDistributionSphericalHeatedMonotonic
     !!{
     Implementation of a heated spherical mass distribution.
     !!}
     private
     class           (massDistributionHeatingClass), pointer     :: massDistributionHeating_ => null()
     double precision                                            :: radiusInitialMinimum               , radiusInitialMaximum, &
          &                                                         radiusFinalMinimum                 , radiusFinalMaximum  , &
          &                                                         radiusVirial
     type            (interpolator                ), allocatable :: massProfile
     logical                                                     :: isBound

   contains
     !![
     <methods>
        <method description="Compute a solution for the heated profile." method="computeSolution" />
     </methods>
     !!]
     final     ::                         sphericalHeatedMonotonicDestructor
     procedure :: computeSolution      => sphericalHeatedMonotonicComputeSolution
     procedure :: density              => sphericalHeatedMonotonicDensity
     procedure :: massEnclosedBySphere => sphericalHeatedMonotonicMassEnclosedBySphere
     procedure :: useUndecorated       => sphericalHeatedMonotonicUseUndecorated
  end type massDistributionSphericalHeatedMonotonic

  interface massDistributionSphericalHeatedMonotonic
     !!{
     Constructors for the \refClass{massDistributionSphericalHeatedMonotonic} mass distribution class.
     !!}
     module procedure sphericalHeatedMonotonicConstructorParameters
     module procedure sphericalHeatedMonotonicConstructorInternal
  end interface massDistributionSphericalHeatedMonotonic

contains

  function sphericalHeatedMonotonicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalHeatedMonotonic} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalHeatedMonotonic)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (massDistributionClass                   ), pointer       :: massDistribution_
    class           (massDistributionHeatingClass            ), pointer       :: massDistributionHeating_
    type            (varying_string                          )                :: nonAnalyticSolver       , componentType, &
         &                                                                       massType
    double precision                                                          :: radiusVirial

    !![
    <inputParameter>
      <name>radiusVirial</name>
      <source>parameters</source>
      <description>The virial radius of the halo.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="massDistribution"        name="massDistribution_"        source="parameters"/>
    <objectBuilder class="massDistributionHeating" name="massDistributionHeating_" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       self=massDistributionSphericalHeatedMonotonic(radiusVirial,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,massDistributionHeating_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"       />
    <objectDestructor name="massDistributionHeating_"/>
    !!]
    return
  end function sphericalHeatedMonotonicConstructorParameters
  
  function sphericalHeatedMonotonicConstructorInternal(radiusVirial,nonAnalyticSolver,massDistribution_,massDistributionHeating_,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalHeatedMonotonic} mass distribution class.
    !!}
    implicit none
    type            (massDistributionSphericalHeatedMonotonic)                          :: self
    double precision                                          , intent(in   )           :: radiusVirial
    class           (massDistributionSpherical               ), intent(in   ), target   :: massDistribution_
    class           (massDistributionHeatingClass            ), intent(in   ), target   :: massDistributionHeating_
    type            (enumerationNonAnalyticSolversType       ), intent(in   )           :: nonAnalyticSolver
    type            (enumerationComponentTypeType            ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType                 ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="radiusVirial, nonAnalyticSolver, *massDistribution_, *massDistributionHeating_, componentType, massType"/>
    !!]
 
    ! Construct the object.
    self%radiusInitialMinimum=+huge(0.0d0)
    self%radiusInitialMaximum=-huge(0.0d0)
    self%radiusFinalMinimum  =+huge(0.0d0)
    self%radiusFinalMaximum  =-huge(0.0d0)
    self%isBound             =.true.
    self%dimensionless       =.false.
    return
  end function sphericalHeatedMonotonicConstructorInternal

  subroutine sphericalHeatedMonotonicDestructor(self)
    !!{
    Destructor for the \refClass{massDistributionSphericalHeatedMonotonic} mass distribution class.
    !!}
    implicit none
    type(massDistributionSphericalHeatedMonotonic), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"       />
    <objectDestructor name="self%massDistributionHeating_"/>
    !!]
    return
  end subroutine sphericalHeatedMonotonicDestructor

  logical function sphericalHeatedMonotonicUseUndecorated(self) result(useUndecorated)
    !!{
    Determines whether to use the undecorated solution.
    !!}
    implicit none
    class(massDistributionSphericalHeatedMonotonic), intent(inout) :: self

    useUndecorated=self%nonAnalyticSolver == nonAnalyticSolversFallThrough .or. self%massDistributionHeating_%specificEnergyIsEverywhereZero()
    return
  end function sphericalHeatedMonotonicUseUndecorated
  
  double precision function sphericalHeatedMonotonicMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (massDistributionSphericalHeatedMonotonic), intent(inout), target :: self
    double precision                                          , intent(in   )         :: radius

    if (self%massDistributionHeating_%specificEnergyIsEverywhereZero()) then
       ! No heating - use the unheated solution.
       mass=self%massDistribution_%massEnclosedBySphere(radius)
    else if (radius <= 0.0d0) then
       ! Non-positive radius, mass must be zero.
       mass=0.0d0
    else
       ! Compute the solution (as needed).
       call self%computeSolution(radius)
       ! For bound halos, interpolate to find the enclosed mass. For unbound halos the enclosed mass is zero.
       if (self%isBound) then
          if (radius < self%radiusFinalMinimum) then
             ! Assume constant density below the minimum radius.
             mass   =+exp(self%massProfile%interpolate(log(self%radiusFinalMinimum))) &
                  &  *(                                                               &
                  &    +     radius                                                   &
                  &    /self%radiusFinalMinimum                                       &
                  &   )**3
          else
             mass=+exp(self%massProfile%interpolate(log(radius)))
          end if
       else
          mass=0.0d0
       end if
    end if
    return
  end function sphericalHeatedMonotonicMassEnclosedBySphere

  double precision function sphericalHeatedMonotonicDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionSphericalHeatedMonotonic), intent(inout) :: self
    class           (coordinate                              ), intent(in   ) :: coordinates
    double precision                                                          :: radius     , radius_
    
    if (self%massDistributionHeating_%specificEnergyIsEverywhereZero()) then
       ! No heating, the density is unchanged.
       density=+self%massDistribution_%density(coordinates)
       return
    end if
    radius=coordinates%rSpherical()
    call self%computeSolution(radius)
    ! For bound halos, interpolate to find the density. For unbound halos the density is zero.
    if (self%isBound) then
       radius_=max(radius,self%radiusFinalMinimum)
       density=+    self%massProfile%derivative (log(radius_))  &
            &  *exp(self%massProfile%interpolate(log(radius_))) &
            &  /4.0d0                                           &
            &  /Pi                                              &
            &  /radius_**3
    else
       density=+0.0d0
    end if
    return
  end function sphericalHeatedMonotonicDensity

  subroutine sphericalHeatedMonotonicComputeSolution(self,radius)
    !!{
    Compute the solution for the heated density profile.
    !!}
    use :: Numerical_Ranges                , only : Make_Range                    , rangeTypeLogarithmic
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Table_Labels                    , only : extrapolationTypeFix          , extrapolationTypeZero
    implicit none
    class           (massDistributionSphericalHeatedMonotonic), intent(inout)             :: self
    double precision                                          , intent(in   )             :: radius
    double precision                                          , parameter                 :: radiusFractionMinimum=1.0d-6, radiusFractionMaximum=10.0d0
    integer                                                   , parameter                 :: countPerDecadeRadius =100
    double precision                                          , allocatable, dimension(:) :: massEnclosed                , massShell                   , &
         &                                                                                   radiusInitial               , radiusFinal                 , &
         &                                                                                   energyFinal                 , perturbation
    logical                                                   , allocatable, dimension(:) :: isBound
    integer                                                                               :: i                           , countRadii

    ! Nothing to do if profile is already tabulated.
    if (.not.self%isBound .or. allocated(self%massProfile)) return
    ! Choose extent of radii at which to tabulate the initial profile.
    self%radiusInitialMinimum=radiusFractionMinimum*self%radiusVirial
    self%radiusInitialMaximum=radiusFractionMaximum*self%radiusVirial
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
       massEnclosed(i)=+self%massDistribution_         %massEnclosedBySphere(radiusInitial(i)                       )
       perturbation(i)=+2.0d0                                                                                         &
            &          *self%massDistributionHeating_  %specificEnergy      (radiusInitial(i),self%massDistribution_) &
            &          /gravitationalConstant_internal                                                                &
            &          /                                                     massEnclosed (i)                         &
            &          *                                                     radiusInitial(i)
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
    energyFinal=+gravitationalConstant_internal &
         &      *massEnclosed                   &
         &      /radiusInitial                  &
         &      *(                              &
         &        -1.0d0                        &
         &        +perturbation                 &
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
       radiusFinal=-gravitationalConstant_internal &
            &      *massEnclosed                   &
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
  end subroutine sphericalHeatedMonotonicComputeSolution
