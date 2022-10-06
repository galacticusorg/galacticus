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

!+    Contributions to this file made by: Xiaolong Du, Andrew Benson.

  !!{
  An implementation of exponentially truncated dark matter halo profiles \cite{kazantzidis_2006}.
  !!}

  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversType, enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOTruncatedExponential">
   <description>exponentially truncated dark matter halo profiles \cite{kazantzidis_2006}.</description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOTruncatedExponential
     !!{
     A dark matter halo profile class implementing exponentially truncated dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_ => null()
     double precision                                             :: radiusFractionalDecay                       , alpha                                                  , &
          &                                                          beta                                        , gamma
     type            (enumerationNonAnalyticSolversType)          :: nonAnalyticSolver
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8                  )           :: lastUniqueID
     ! Stored values of computed quantities.
     double precision                                             :: radialVelocityDispersionVirialRadiusPrevious, radialVelocityDispersionVirialRadiusUntruncatedPrevious, &
          &                                                          enclosingMassRadiusPrevious                 , kappaPrevious                                          , &
          &                                                          gammaFunctionIncompletePrevious             , densityNormalizationPrevious                           , &
          &                                                          massVirialPrevious
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
       <method description="Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in units of Mpc)." method="truncationFunction" />
     </methods>
     !!]
     final     ::                                      truncatedExponentialDestructor
     procedure :: autoHook                          => truncatedExponentialAutoHook
     procedure :: calculationReset                  => truncatedExponentialCalculationReset
     procedure :: density                           => truncatedExponentialDensity
     procedure :: densityLogSlope                   => truncatedExponentialDensityLogSlope
     procedure :: radiusEnclosingDensity            => truncatedExponentialRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => truncatedExponentialRadiusEnclosingMass
     procedure :: radialMoment                      => truncatedExponentialRadialMoment
     procedure :: enclosedMass                      => truncatedExponentialEnclosedMass
     procedure :: potential                         => truncatedExponentialPotential
     procedure :: circularVelocity                  => truncatedExponentialCircularVelocity
     procedure :: radiusCircularVelocityMaximum     => truncatedExponentialRadiusCircularVelocityMaximum
     procedure :: circularVelocityMaximum           => truncatedExponentialCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => truncatedExponentialRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => truncatedExponentialRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => truncatedExponentialRotationNormalization
     procedure :: energy                            => truncatedExponentialEnergy
     procedure :: kSpace                            => truncatedExponentialKSpace
     procedure :: freefallRadius                    => truncatedExponentialFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => truncatedExponentialFreefallRadiusIncreaseRate
     procedure :: truncationFunction                => truncatedExponentialTruncationFunction
  end type darkMatterProfileDMOTruncatedExponential

  interface darkMatterProfileDMOTruncatedExponential
     !!{
     Constructors for the {\normalfont \ttfamily exponentially truncated} dark matter halo profile class.
     !!}
     module procedure truncatedExponentialConstructorParameters
     module procedure truncatedExponentialConstructorInternal
  end interface darkMatterProfileDMOTruncatedExponential

contains

  function truncatedExponentialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily exponentially truncated} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOTruncatedExponential)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass               ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    type            (varying_string                          )                :: nonAnalyticSolver
    double precision                                                          :: radiusFractionalDecay, alpha, &
         &                                                                       beta                 , gamma

    !![
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusFractionalDecay</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <description>The truncation scale (in units of the virial radius).</description>
    </inputParameter>
    <inputParameter>
      <name>alpha</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <description>Parameter $\alpha$ in the \cite{kazantzidis_2006} truncated profile.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <defaultValue>3.0d0</defaultValue>
      <source>parameters</source>
      <description>Parameter $\beta$ in the \cite{kazantzidis_2006} truncated profile.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <description>Parameter $\gamma$ in the \cite{kazantzidis_2006} truncated profile.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=darkMatterProfileDMOTruncatedExponential(radiusFractionalDecay,alpha,beta,gamma,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function truncatedExponentialConstructorParameters

  function truncatedExponentialConstructorInternal(radiusFractionalDecay,alpha,beta,gamma,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily exponentially truncated} dark matter profile class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (darkMatterProfileDMOTruncatedExponential)                        :: self
    class           (darkMatterProfileDMOClass               ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                          , intent(in   )         :: radiusFractionalDecay, alpha, &
         &                                                                               beta                 , gamma
    type            (enumerationNonAnalyticSolversType       ), intent(in   )         :: nonAnalyticSolver
    !![
    <constructorAssign variables="radiusFractionalDecay,alpha,beta,gamma,nonAnalyticSolver,*darkMatterProfileDMO_,*darkMatterHaloScale_"/>
    !!]

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    self%lastUniqueID                                           =-1_kind_int8
    self%genericLastUniqueID                                    =-1_kind_int8
    self%kappaPrevious                                          =-huge(0.0d0)
    self%enclosingMassRadiusPrevious                            =-1.0d0
    self%radialVelocityDispersionVirialRadiusPrevious           =-1.0d0
    self%radialVelocityDispersionVirialRadiusUntruncatedPrevious=-1.0d0
    return
  end function truncatedExponentialConstructorInternal

  subroutine truncatedExponentialAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOTruncatedExponential), intent(inout) :: self

    call calculationResetEvent%attach(self,truncatedExponentialCalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileDMOTruncatedExponential')
    return
  end subroutine truncatedExponentialAutoHook

  subroutine truncatedExponentialDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily exponentially truncated} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMOTruncatedExponential), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    if (calculationResetEvent%isAttached(self,truncatedExponentialCalculationReset)) call calculationResetEvent%detach(self,truncatedExponentialCalculationReset)
    return
  end subroutine truncatedExponentialDestructor

  subroutine truncatedExponentialCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node

    self%lastUniqueID                                           =node%uniqueID()
    self%genericLastUniqueID                                    =node%uniqueID()
    self%kappaPrevious                                          =-huge(0.0d0)
    self%enclosingMassRadiusPrevious                            =-1.0d0
    self%radialVelocityDispersionVirialRadiusPrevious           =-1.0d0
    self%radialVelocityDispersionVirialRadiusUntruncatedPrevious=-1.0d0
    self%genericEnclosedMassRadiusMinimum                       =+huge(0.0d0)
    self%genericEnclosedMassRadiusMaximum                       =-huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMinimum           =+huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMaximum           =-huge(0.0d0)
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine truncatedExponentialCalculationReset

  subroutine truncatedExponentialTruncationFunction(self,node,radius,multiplier,multiplierGradient)
    !!{
    Return the scaled truncation radial coordinate, and the truncation multiplier.
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout)           :: self
    type            (treeNode                                ), intent(inout)           :: node
    double precision                                          , intent(in   )           :: radius
    double precision                                          , intent(  out), optional :: multiplier  , multiplierGradient
    double precision                                                                    :: radiusVirial, radiusDecay       , &
         &                                                                                 multiplier_

    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    radiusVirial=self%darkMatterHaloScale_%radiusVirial(node)
    if (radius <= radiusVirial) then
       if (present(multiplier        )) multiplier        =+1.0d0
       if (present(multiplierGradient)) multiplierGradient=+0.0d0
    else
       if (self%kappaPrevious == -huge(0.0d0)) call recomputeKappa(self,node)
       radiusDecay                                        =+self%radiusFractionalDecay*radiusVirial
       multiplier_                                        =+self%darkMatterProfileDMO_%density(node,radiusVirial) &
            &                                              /self%darkMatterProfileDMO_%density(node,radius      ) &
            &                                              *(                                                     &
            &                                                +radius                                              &
            &                                                /radiusVirial                                        &
            &                                               )**self%kappaPrevious                                 &
            &                                              *exp(                                                  &
            &                                                   -(                                                &
            &                                                     +radius                                         &
            &                                                     -radiusVirial                                   &
            &                                                    )                                                &
            &                                                   /  radiusDecay                                    &
            &                                                  )
       if (present(multiplier        )) multiplier        =+multiplier_
       if (present(multiplierGradient)) multiplierGradient=+multiplier_                                                           &
            &                                              *(                                                                     &
            &                                                +self%kappaPrevious                                     /radius      &
            &                                                -1.0d0                                                  /radiusDecay &
            &                                                -self%darkMatterProfileDMO_%densityLogSlope(node,radius)/radius      &
            &                                               )
    end if
    return
  end subroutine truncatedExponentialTruncationFunction

  subroutine recomputeKappa (self,node)
    !!{
    Recompute parameter kappa in the truncation funciton.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile        , treeNode
    use :: Gamma_Functions , only : Gamma_Function_Incomplete_Unnormalized
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    class           (nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile
    double precision                                                          :: radiusVirial     , scaleRadius, &
         &                                                                       concentration

    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    radiusVirial                        =  self             %darkMatterHaloScale_%radiusVirial(node             )
    darkMatterProfile                   => node             %darkMatterProfile                (autoCreate=.true.)
    scaleRadius                         =  darkMatterProfile%scale                            (                 )
    concentration                       =  +radiusVirial                &
         &                                 /scaleRadius
    self%kappaPrevious                  =  -(                           &
         &                                   +self%gamma                &
         &                                   +self%beta                 &
         &                                   *concentration**self%alpha &
         &                                  )                           &
         &                                 /(                           &
         &                                   +1.0d0                     &
         &                                   +concentration**self%alpha &
         &                                  )                           &
         &                                 +1.0d0                       &
         &                                 /self%radiusFractionalDecay
    self%massVirialPrevious             =+self%darkMatterProfileDMO_%enclosedMass(node                    ,radiusVirial                    )
    self%densityNormalizationPrevious   =+4.0d0                                                                                              &
         &                               *Pi                                                                                                 &
         &                               *self%darkMatterProfileDMO_%density     (node                    ,radiusVirial                    ) &
         &                               *radiusVirial**3                                                                                    &
         &                               *self%radiusFractionalDecay**           (3.0d0+self%kappaPrevious                                 ) &
         &                               *exp                                    (                         1.0d0/self%radiusFractionalDecay)
    self%gammaFunctionIncompletePrevious=+Gamma_Function_Incomplete_Unnormalized (3.0d0+self%kappaPrevious,1.0d0/self%radiusFractionalDecay)
    return
  end subroutine recomputeKappa

  double precision function truncatedExponentialDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: radius
    double precision                                                          :: multiplier

    call self%truncationFunction(node,radius,multiplier=multiplier)
    truncatedExponentialDensity=+self%darkMatterProfileDMO_%density(node,radius) &
         &                      *multiplier
    return
  end function truncatedExponentialDensity

  double precision function truncatedExponentialDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: radius
    double precision                                                          :: multiplier, multiplierGradient

    call self%truncationFunction(node,radius,multiplier=multiplier,multiplierGradient=multiplierGradient)
    if (multiplier > 0.0d0) then
       truncatedExponentialDensityLogSlope=+self%darkMatterProfileDMO_%densityLogSlope(node,radius) &
            &                              +radius                                                  &
            &                              *multiplierGradient                                      &
            &                              /multiplier
    else
       truncatedExponentialDensityLogSlope=+0.0d0
    end if
    return
  end function truncatedExponentialDensityLogSlope

  double precision function truncatedExponentialRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout), target :: self
    type            (treeNode                                ), intent(inout), target :: node
    double precision                                          , intent(in   )         :: density

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity         (node,density)
    else
       truncatedExponentialRadiusEnclosingDensity=self                      %radiusEnclosingDensityNumerical(node,density)
    end if
    return
  end function truncatedExponentialRadiusEnclosingDensity

  double precision function truncatedExponentialRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    double precision                                          , intent(in   )          :: mass
    class           (nodeComponentBasic                      ),                pointer :: basic
    double precision                                                                   :: massVirial

    basic      => node %basic()
    massVirial =  basic%mass ()
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough .or. mass <= massVirial) then
       truncatedExponentialRadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass         (node,mass)
    else
       truncatedExponentialRadiusEnclosingMass=self                      %radiusEnclosingMassNumerical(node,mass)
    end if
    return
  end function truncatedExponentialRadiusEnclosingMass

  double precision function truncatedExponentialRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout)           :: self
    type            (treeNode                                ), intent(inout)           :: node
    double precision                                          , intent(in   )           :: moment
    double precision                                          , intent(in   ), optional :: radiusMinimum, radiusMaximum

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialRadialMoment=self%darkMatterProfileDMO_%radialMoment         (node,moment,radiusMinimum,radiusMaximum)
    else
       truncatedExponentialRadialMoment=self                      %radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    end if
    return
  end function truncatedExponentialRadialMoment

  double precision function truncatedExponentialEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Gamma_Functions         , only : Gamma_Function_Incomplete_Unnormalized
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: radius
    double precision                                                          :: radiusVirial, radiusDecay

    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    radiusVirial=self%darkMatterHaloScale_%radiusVirial(node)
    if (radius <= radiusVirial) then
       truncatedExponentialEnclosedMass=+self%darkMatterProfileDMO_%enclosedMass(node,radius      )
    else
       if (self%kappaPrevious == -huge(0.0d0)) call recomputeKappa(self,node)
       radiusDecay                     =+self%radiusFractionalDecay*radiusVirial
       truncatedExponentialEnclosedMass=+self%massVirialPrevious                                                               &
            &                           +self%densityNormalizationPrevious                                                     &
            &                           *(                                                                                     &
            &                             +self%gammaFunctionIncompletePrevious                                                &
            &                             -Gamma_Function_Incomplete_Unnormalized(3.0d0+self%kappaPrevious,radius/radiusDecay) &
            &                            )
    end if
    return
  end function truncatedExponentialEnclosedMass

  double precision function truncatedExponentialPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout)           :: self
    type            (treeNode                                ), intent(inout), target   :: node
    double precision                                          , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType       ), intent(  out), optional :: status

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialPotential=self%darkMatterProfileDMO_%potential         (node,radius,status)
    else
       truncatedExponentialPotential=self                      %potentialNumerical(node,radius,status)
    end if
    return
  end function truncatedExponentialPotential

  double precision function truncatedExponentialCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: radius

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialCircularVelocity=self%darkMatterProfileDMO_%circularVelocity         (node,radius)
    else
       truncatedExponentialCircularVelocity=self                      %circularVelocityNumerical(node,radius)
    end if
    return
  end function truncatedExponentialCircularVelocity

  double precision function truncatedExponentialRadiusCircularVelocityMaximum(self,node)
    !!{
    Returns the radius (in Mpc) at which the maximum circular velocity is acheived in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialRadiusCircularVelocityMaximum=self%darkMatterProfileDMO_%radiusCircularVelocityMaximum         (node)
    else
       truncatedExponentialRadiusCircularVelocityMaximum=self                      %radiusCircularVelocityMaximumNumerical(node)
    end if
    return
  end function truncatedExponentialRadiusCircularVelocityMaximum

  double precision function truncatedExponentialCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialCircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum         (node)
    else
       truncatedExponentialCircularVelocityMaximum=self                      %circularVelocityMaximumNumerical(node)
    end if
    return
  end function truncatedExponentialCircularVelocityMaximum

  double precision function truncatedExponentialRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: radius
    double precision                                                          :: radiusVirial

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialRadialVelocityDispersion=self%darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
    else
       if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
       radiusVirial=self%darkMatterHaloScale_%radiusVirial(node)
       if (radius >= radiusVirial) then
          truncatedExponentialRadialVelocityDispersion=self%radialVelocityDispersionNumerical(node,radius)
       else
          if (self%radialVelocityDispersionVirialRadiusPrevious < 0.0d0 .or. self%radialVelocityDispersionVirialRadiusUntruncatedPrevious < 0.0d0) then
             self%radialVelocityDispersionVirialRadiusPrevious           =self%radialVelocityDispersionNumerical             (node,radiusVirial)
             self%radialVelocityDispersionVirialRadiusUntruncatedPrevious=self%darkMatterProfileDMO_%radialVelocityDispersion(node,radiusVirial)
          end if
          truncatedExponentialRadialVelocityDispersion=sqrt(                                                                             &
               &                                            +self%darkMatterProfileDMO_%radialVelocityDispersion  (node,radius      )**2 &
               &                                            +self%density                                         (node,radiusVirial)    &
               &                                            /self%density                                         (node,radius      )    &
               &                                            *(                                                                           &
               &                                              +self%radialVelocityDispersionVirialRadiusPrevious                     **2 &
               &                                              -self%radialVelocityDispersionVirialRadiusUntruncatedPrevious          **2 &
               &                                             )                                                                           &
               &                                           )
       end if
    end if
    return
  end function truncatedExponentialRadialVelocityDispersion

  double precision function truncatedExponentialRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: specificAngularMomentum

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum         (node,specificAngularMomentum)
    else
       truncatedExponentialRadiusFromSpecificAngularMomentum=self                      %radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    end if
    return
  end function truncatedExponentialRadiusFromSpecificAngularMomentum

  double precision function truncatedExponentialRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class(darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialRotationNormalization=self%darkMatterProfileDMO_%rotationNormalization         (node)
    else
       truncatedExponentialRotationNormalization=self                      %rotationNormalizationNumerical(node)
    end if
    return
  end function truncatedExponentialRotationNormalization

  double precision function truncatedExponentialEnergy(self,node)
    !!{
    Return the energy of a truncatedExponential halo density profile.
    !!}
    implicit none
    class(darkMatterProfileDMOTruncatedExponential), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialEnergy=self%darkMatterProfileDMO_%energy         (node)
    else
       truncatedExponentialEnergy=self                      %energyNumerical(node)
    end if
    return
  end function truncatedExponentialEnergy

  double precision function truncatedExponentialKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the truncatedExponential density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout)         :: self
    type            (treeNode                                ), intent(inout), target :: node
    double precision                                          , intent(in   )         :: waveNumber

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialKSpace=self%darkMatterProfileDMO_%kSpace         (node,waveNumber)
    else
       truncatedExponentialKSpace=self                      %kSpaceNumerical(node,waveNumber)
    end if
    return
  end function truncatedExponentialKSpace

  double precision function truncatedExponentialFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the truncatedExponential density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout), target :: self
    type            (treeNode                                ), intent(inout), target :: node
    double precision                                          , intent(in   )         :: time

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialFreefallRadius=self%darkMatterProfileDMO_%freefallRadius         (node,time)
    else
       truncatedExponentialFreefallRadius=self                      %freefallRadiusNumerical(node,time)
    end if
    return
  end function truncatedExponentialFreefallRadius

  double precision function truncatedExponentialFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the truncatedExponential density profile at the specified {\normalfont
    \ttfamily time} (given in Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout), target :: self
    type            (treeNode                                ), intent(inout), target :: node
    double precision                                          , intent(in   )         :: time

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       truncatedExponentialFreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate         (node,time)
    else
       truncatedExponentialFreefallRadiusIncreaseRate=self                      %freefallRadiusIncreaseRateNumerical(node,time)
    end if
    return
  end function truncatedExponentialFreefallRadiusIncreaseRate
