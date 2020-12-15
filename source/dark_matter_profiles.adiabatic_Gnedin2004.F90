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

  !% An implementation of adiabaticGnedin2004 dark matter halo profiles.

  use :: Cosmology_Parameters        , only : cosmologyParameters                , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScale                , darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO    , only : darkMatterProfileDMO               , darkMatterProfileDMOClass
  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough
  use :: Galactic_Structure_Options  , only : componentTypeAll                   , massTypeBaryonic                    , radiusLarge                  , weightByMass, &
          &                                   weightIndexNull
  use :: Math_Exponentiation         , only : fastExponentiator

  ! Number of previous radius solutions to store.
  !# <scoping>
  !#  <module variables="adiabaticGnedin2004StoreCount"/>
  !# </scoping>
  integer, parameter :: adiabaticGnedin2004StoreCount=10

  !# <darkMatterProfile name="darkMatterProfileAdiabaticGnedin2004">
  !#  <description>
  !#   A dark matter profile class which applies adiabatic contraction to the halo as it responds to the presence of
  !#   baryons. Adiabatic contraction follows the algorithm of \cite{gnedin_response_2004}. The parameters $A$ and $\omega$ of
  !#   that model are specified via input parameters {\normalfont \ttfamily A} and {\normalfont \ttfamily omega} respectively.
  !#
  !#   Given the final radius, $r_\mathrm{f}$, the corresponding initial radius, $r_\mathrm{i}$, is found by solving:
  !#   \begin{equation}
  !#   f_\mathrm{i} M_\mathrm{total,0}(\bar{r}_\mathrm{i}) r_\mathrm{i} = f_\mathrm{f} M_\mathrm{total,0}(\bar{r}_\mathrm{i})
  !#   r_\mathrm{f} + V^2_\mathrm{b}(\bar{r}_\mathrm{f}) \bar{r}_\mathrm{f} r_\mathrm{f}/ \mathrm{G},
  !#    \label{eq:adiabaticContractionGnedinSolution}
  !#   \end{equation}
  !#   where $M_\mathrm{total,0}(r)$ is the initial total matter profile, $V_\mathrm{b}(r)$ is the baryonic contribution to the
  !#   rotation curve, $f_\mathrm{i}$, is the fraction of mass within the virial radius compared to the node mass\footnote{In
  !#   \protect\glc\ the ``node mass'' refers to the total mass of the node, assuming it has the universal complement of
  !#   baryons. Since some halos may contain less than the complete complement of baryons it is possible that $f_\mathrm{i}&lt;1$.},
  !#   $f_\mathrm{f}=(\Omega_\mathrm{M}-\Omega_\mathrm{b})/\Omega_\mathrm{M}+M_\mathrm{satellite, baryonic}/M_\mathrm{total}$,
  !#   $M_\mathrm{satellite, baryonic}$ is the baryonic mass in any satellite halos, $M_\mathrm{total}$ is the node mass, and
  !#   \begin{equation}
  !#   {\bar{r} \over r_\mathrm{vir}} = A \left({r \over r_\mathrm{vir}}\right)^\omega,
  !#   \end{equation}
  !#   where $r_\mathrm{vir}$ is the virial radius. Note that we explicitly assume that the initial, uncontracted total density
  !#   profile has the same shape as the initial dark matter density profile, that contraction of the halo occurs with no shell
  !#   crossing, and that satellite halos trace the dark matter profile of their host halo.
  !#   The derivative, $\mathrm{d} r_\mathrm{f}/\mathrm{d}d_\mathrm{i}\equiv r^\prime_\mathrm{i}$ is found by taking the
  !#   derivative of eqn.~(\ref{eq:adiabaticContractionGnedinSolution}) to give:
  !#   \begin{eqnarray}
  !#    &amp; &amp; f_\mathrm{i} M_\mathrm{total,0}(\bar{r}_\mathrm{i}) r^\prime_\mathrm{i} + f_\mathrm{i} 4 \pi
  !#    \bar{r}_\mathrm{i}^2 \rho_\mathrm{total,0}(\bar{r}_\mathrm{i}) {\mathrm{d} \bar{r}_\mathrm{i}\over\mathrm{d} r_\mathrm{i}}
  !#    r_\mathrm{i} r^\prime_\mathrm{i} \nonumber \\
  !#    &amp; = &amp; f_\mathrm{f} M_\mathrm{total,0}(\bar{r}_\mathrm{i}) + f_\mathrm{i} 4 \pi \bar{r}_\mathrm{i}^2
  !#    \rho_\mathrm{total,0}(\bar{r}_\mathrm{i}) {\mathrm{d} \bar{r}_\mathrm{i}\over\mathrm{d} r_\mathrm{i}} r_\mathrm{f}
  !#    r^\prime_\mathrm{i} \nonumber \\
  !#    &amp; + &amp; V^2_\mathrm{b}(\bar{r}_\mathrm{f}) \bar{r}_\mathrm{f} / \mathrm{G} + V^2_\mathrm{b}(\bar{r}_\mathrm{f})
  !#    {\mathrm{d}\bar{r}_\mathrm{f}\over \mathrm{d} r_\mathrm{f}} r_\mathrm{f}/ \mathrm{G} +
  !#    {\mathrm{d}V^2_\mathrm{b}\over\mathrm{d} \bar{r}_\mathrm{f}}(\bar{r}_\mathrm{f}) {\mathrm{d}\bar{r}_\mathrm{f}\over
  !#    \mathrm{d} r_\mathrm{f}} \bar{r}_\mathrm{f} r_\mathrm{f}/ \mathrm{G},
  !#   \end{eqnarray}
  !#   where
  !#   \begin{equation}
  !#    {\mathrm{d}\bar{r} \over \mathrm{d} r} = A \left({r \over r_\mathrm{vir}}\right)^{\omega-1},
  !#   \end{equation}
  !#   and which can then be solved numerically for $r^\prime_\mathrm{i}$.
  !#  </description>
  !# </darkMatterProfile>
  type, extends(darkMatterProfileClass) :: darkMatterProfileAdiabaticGnedin2004
     !% A dark matter halo profile class implementing adiabaticGnedin2004 dark matter halos.
     private
     class           (cosmologyParametersClass ), pointer                                  :: cosmologyParameters_            => null()
     class           (darkMatterProfileDMOClass), pointer                                  :: darkMatterProfileDMO_           => null()
     integer                                                                               :: nonAnalyticSolver
     ! Parameters of the adiabatic contraction algorithm.
     double precision                                                                      :: A                                        , omega
     ! Stored solutions for reuse.
     integer         (kind=kind_int8           )                                           :: lastUniqueID
     integer                                                                               :: radiusPreviousIndex                      , radiusPreviousIndexMaximum
     double precision                           , dimension(adiabaticGnedin2004StoreCount) :: radiusPrevious                           , radiusInitialPrevious
     type            (fastExponentiator        )                                           :: radiusExponentiator
     ! Quantities used in solving the initial radius root function.
     double precision                                                                      :: baryonicFinalTerm                        , baryonicFinalTermDerivative, &
          &                                                                                   darkMatterDistributedFraction            , initialMassFraction        , &
          &                                                                                   radiusFinal                              , radiusFinalMean            , &
          &                                                                                   radiusFinalMeanSelfDerivative            , radiusInitial_             , &
          &                                                                                   radiusInitialMeanSelfDerivative          , radiusShared               , &
          &                                                                                   radiusVirial                             , darkMatterFraction
     logical                                                                               :: massesComputed
   contains
     !# <methods>
     !#   <method description="Reset memoized calculations." method="calculationReset" />
     !#   <method description="Compute factors needed for solving adiabatic contraction." method="computeFactors" />
     !#   <method description="Compute the orbit-averaged radius for dark matter." method="radiusOrbitalMean" />
     !#   <method description="Compute the derivative of the orbit-averaged radius for dark matter." method="radiusOrbitalMeanDerivative" />
     !#   <method description="Compute the initial radius in the dark matter profile." method="radiusInitial" />
     !#   <method description="Compute the derivative of the initial radius in the dark matter profile." method="radiusInitialDerivative" />
     !# </methods>
     final     ::                                      adiabaticGnedin2004Destructor
     procedure :: autoHook                          => adiabaticGnedin2004AutoHook
     procedure :: calculationReset                  => adiabaticGnedin2004CalculationReset
     procedure :: density                           => adiabaticGnedin2004Density
     procedure :: densityLogSlope                   => adiabaticGnedin2004DensityLogSlope
     procedure :: radiusEnclosingDensity            => adiabaticGnedin2004RadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => adiabaticGnedin2004RadiusEnclosingMass
     procedure :: radialMoment                      => adiabaticGnedin2004RadialMoment
     procedure :: enclosedMass                      => adiabaticGnedin2004EnclosedMass
     procedure :: potential                         => adiabaticGnedin2004Potential
     procedure :: circularVelocity                  => adiabaticGnedin2004CircularVelocity
     procedure :: circularVelocityMaximum           => adiabaticGnedin2004CircularVelocityMaximum
     procedure :: radialVelocityDispersion          => adiabaticGnedin2004RadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => adiabaticGnedin2004RadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => adiabaticGnedin2004RotationNormalization
     procedure :: energy                            => adiabaticGnedin2004Energy
     procedure :: energyGrowthRate                  => adiabaticGnedin2004EnergyGrowthRate
     procedure :: kSpace                            => adiabaticGnedin2004KSpace
     procedure :: freefallRadius                    => adiabaticGnedin2004FreefallRadius
     procedure :: freefallRadiusIncreaseRate        => adiabaticGnedin2004FreefallRadiusIncreaseRate
     procedure :: radiusInitial                     => adiabaticGnedin2004RadiusInitial
     procedure :: radiusInitialDerivative           => adiabaticGnedin2004RadiusInitialDerivative
     procedure :: computeFactors                    => adiabaticGnedin2004ComputeFactors
     procedure :: radiusOrbitalMean                 => adiabaticGnedin2004RadiusOrbitalMean
     procedure :: radiusOrbitalMeanDerivative       => adiabaticGnedin2004RadiusOrbitalMeanDerivative
  end type darkMatterProfileAdiabaticGnedin2004

  interface darkMatterProfileAdiabaticGnedin2004
     !% Constructors for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter halo profile class.
     module procedure adiabaticGnedin2004ConstructorParameters
     module procedure adiabaticGnedin2004ConstructorInternal
  end interface darkMatterProfileAdiabaticGnedin2004

  ! Module-scope quantities used in solving the initial radius root function.
  integer                                               , parameter :: adiabaticGnedin2004ComponentType=componentTypeAll, adiabaticGnedin2004MassType   =massTypeBaryonic, &
       &                                                               adiabaticGnedin2004WeightBy     =weightByMass    , adiabaticGnedin2004WeightIndex=weightIndexNull
  type            (treeNode                            ), pointer   :: adiabaticGnedin2004Node
  class           (darkMatterProfileAdiabaticGnedin2004), pointer   :: adiabaticGnedin2004Self
  !$omp threadprivate(adiabaticGnedin2004Self,adiabaticGnedin2004Node)

contains

  function adiabaticGnedin2004ConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter halo profile class.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileAdiabaticGnedin2004)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (cosmologyParametersClass            ), pointer       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass            ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass           ), pointer       :: darkMatterProfileDMO_
    type            (varying_string                      )                :: nonAnalyticSolver
    double precision                                                      :: A                    , omega

    !# <inputParameter>
    !#   <name>A</name>
    !#   <defaultSource>(\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultSource>
    !#   <defaultValue>0.80d0</defaultValue>
    !#   <description>The parameter $A$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>omega</name>
    !#   <defaultSource>(\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultSource>
    !#   <defaultValue>0.77d0</defaultValue>
    !#   <description>The parameter $\omega$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>nonAnalyticSolver</name>
    !#   <defaultValue>var_str('fallThrough')</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring adiabatic contraction by baryons is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !# <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    self=darkMatterProfileAdiabaticGnedin2004(A,omega,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_" />
    !# <objectDestructor name="darkMatterHaloScale_" />
    !# <objectDestructor name="darkMatterProfileDMO_"/>
    return
  end function adiabaticGnedin2004ConstructorParameters

  function adiabaticGnedin2004ConstructorInternal(A,omega,nonAnalyticSolver,cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !% Generic constructor for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter profile class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (darkMatterProfileAdiabaticGnedin2004)                        :: self
    double precision                                      , intent(in   )         :: A                   , omega
    class           (cosmologyParametersClass            ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterProfileDMOClass           ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), intent(in   ), target :: darkMatterHaloScale_
    integer                                               , intent(in   )         :: nonAnalyticSolver
    !# <constructorAssign variables="A, omega, nonAnalyticSolver, *cosmologyParameters_, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Galacticus_Error_Report('invalid non-analytic solver type'//{introspection:location})
    ! Construct the object.
    self%lastUniqueID       =-1_kind_int8
    self%massesComputed     =.false.
    self%radiusExponentiator=fastExponentiator(1.0d-3,1.0d0,omega,1.0d4,.false.)
    ! Evaluate the dark matter fraction.
    self%darkMatterFraction=+1.0d0                                   &
         &                  -self%cosmologyParameters_%OmegaBaryon() &
         &                  /self%cosmologyParameters_%OmegaMatter()
    return
  end function adiabaticGnedin2004ConstructorInternal

  subroutine adiabaticGnedin2004AutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self

    call calculationResetEvent%attach(self,adiabaticGnedin2004CalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine adiabaticGnedin2004AutoHook

  subroutine adiabaticGnedin2004Destructor(self)
    !% Destructor for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter halo profile class.
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_" />
    !# <objectDestructor name="self%darkMatterHaloScale_" />
    !# <objectDestructor name="self%darkMatterProfileDMO_"/>
    call calculationResetEvent%detach(self,adiabaticGnedin2004CalculationReset)
    return
  end subroutine adiabaticGnedin2004Destructor

  subroutine adiabaticGnedin2004CalculationReset(self,node)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    ! Reset calculations for this profile.
    self%lastUniqueID              =node%uniqueID()
    self%radiusPreviousIndex       = 0
    self%radiusPreviousIndexMaximum= 0
    self%radiusPrevious            =-1.0d0
    self%massesComputed            =.false.
    return
  end subroutine adiabaticGnedin2004CalculationReset

  double precision function adiabaticGnedin2004EnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius

    adiabaticGnedin2004EnclosedMass=+self%darkMatterFraction                                                       &
         &                          *self%darkMatterProfileDMO_%enclosedMass(node,self%radiusInitial(node,radius))
    return
  end function adiabaticGnedin2004EnclosedMass

  double precision function adiabaticGnedin2004Density(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    double precision                                                      :: radiusInitial , radiusInitialDerivative, &
         &                                                                   densityInitial

    radiusInitial =self                      %radiusInitial(node,radius       )
    densityInitial=self%darkMatterProfileDMO_%density      (node,radiusInitial)
    if (radius == radiusInitial) then
       adiabaticGnedin2004Density=+self%darkMatterFraction &
            &                     *densityInitial
    else
       radiusInitialDerivative   = self%radiusInitialDerivative(node,radius)
       adiabaticGnedin2004Density=+self%darkMatterFraction    &
            &                     *densityInitial             &
            &                     *(                          &
            &                       +radiusInitial            &
            &                       /radius                   &
            &                      )                      **2 &
            &                     *radiusInitialDerivative
    end if
    return
  end function adiabaticGnedin2004Density

  double precision function adiabaticGnedin2004DensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    double precision                                                      :: radiusInitial

    radiusInitial=self%radiusInitial(node,radius)
    if (radius == radiusInitial .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004DensityLogSlope=self%darkMatterProfileDMO_%densityLogSlope         (node,radius)
    else
       adiabaticGnedin2004DensityLogSlope=self                      %densityLogSlopeNumerical(node,radius)
    end if
    return
  end function adiabaticGnedin2004DensityLogSlope

  double precision function adiabaticGnedin2004RadiusEnclosingDensity(self,node,density)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: density

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity         (node,density)
    else
       adiabaticGnedin2004RadiusEnclosingDensity=self                      %radiusEnclosingDensityNumerical(node,density)
    end if
    return
  end function adiabaticGnedin2004RadiusEnclosingDensity

  double precision function adiabaticGnedin2004RadiusEnclosingMass(self,node,mass)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: mass

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass         (node,mass)
    else
       adiabaticGnedin2004RadiusEnclosingMass=self                      %radiusEnclosingMassNumerical(node,mass)
    end if
    return
  end function adiabaticGnedin2004RadiusEnclosingMass

  double precision function adiabaticGnedin2004RadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the radial moment of the density in the dark matter profile of {\normalfont \ttfamily node} between the given
    !% {\normalfont \ttfamily radiusMinimum} and {\normalfont \ttfamily radiusMaximum} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    double precision                                      , intent(in   )           :: moment
    double precision                                      , intent(in   ), optional :: radiusMinimum, radiusMaximum

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadialMoment=self%darkMatterProfileDMO_%radialMoment         (node,moment,radiusMinimum,radiusMaximum)
    else
       adiabaticGnedin2004RadialMoment=self                      %radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    end if
    return
  end function adiabaticGnedin2004RadialMoment

  double precision function adiabaticGnedin2004Potential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target   :: self
    type            (treeNode                            ), intent(inout), target   :: node
    double precision                                      , intent(in   )           :: radius
    integer                                               , intent(  out), optional :: status
    double precision                                                                :: radiusInitial

    radiusInitial=self%radiusInitial(node,radius)
    if (radius == radiusInitial .or. self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       ! No adiabatic contraction - use the dark-matter-only result.
       adiabaticGnedin2004Potential=+self%darkMatterFraction                                           &
            &                       *self%darkMatterProfileDMO_%potential         (node,radius,status)
    else
       ! Adiabatic contraction is present - fall back to using a numerical calculation.
       adiabaticGnedin2004Potential=+self                      %potentialNumerical(node,radius,status)
    end if
    return
  end function adiabaticGnedin2004Potential

  double precision function adiabaticGnedin2004CircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004CircularVelocity=self%darkMatterProfileDMO_%circularVelocity         (node,radius)
    else
       adiabaticGnedin2004CircularVelocity=self                      %circularVelocityNumerical(node,radius)
    end if
    return
  end function adiabaticGnedin2004CircularVelocity

  double precision function adiabaticGnedin2004CircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004CircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum         (node)
    else
       adiabaticGnedin2004CircularVelocityMaximum=self                      %circularVelocityMaximumNumerical(node)
    end if
    return
  end function adiabaticGnedin2004CircularVelocityMaximum

  double precision function adiabaticGnedin2004RadialVelocityDispersion(self,node,radius)
    !% Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadialVelocityDispersion=self%darkMatterProfileDMO_%radialVelocityDispersion         (node,radius)
    else
       adiabaticGnedin2004RadialVelocityDispersion=self                      %radialVelocityDispersionNumerical(node,radius)
    end if
    return
  end function adiabaticGnedin2004RadialVelocityDispersion

  double precision function adiabaticGnedin2004RadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target :: self
    type            (treeNode                            ), intent(inout)         :: node
    double precision                                      , intent(in   )         :: specificAngularMomentum

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum         (node,specificAngularMomentum)
    else
       adiabaticGnedin2004RadiusFromSpecificAngularMomentum=self                      %radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    end if
    return
  end function adiabaticGnedin2004RadiusFromSpecificAngularMomentum

  double precision function adiabaticGnedin2004RotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004RotationNormalization=self%darkMatterProfileDMO_%rotationNormalization         (node)
    else
       adiabaticGnedin2004RotationNormalization=self                      %rotationNormalizationNumerical(node)
    end if
    return
  end function adiabaticGnedin2004RotationNormalization

  double precision function adiabaticGnedin2004Energy(self,node)
    !% Return the energy of a adiabaticGnedin2004 halo density profile.
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004Energy=self%darkMatterProfileDMO_%energy         (node)
    else
       adiabaticGnedin2004Energy=self                      %energyNumerical(node)
    end if
    return
  end function adiabaticGnedin2004Energy

  double precision function adiabaticGnedin2004EnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of a adiabaticGnedin2004 halo density profile.
    implicit none
    class(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004EnergyGrowthRate=self%darkMatterProfileDMO_%energyGrowthRate         (node)
    else
       adiabaticGnedin2004EnergyGrowthRate=self                      %energyGrowthRateNumerical(node)
    end if
    return
  end function adiabaticGnedin2004EnergyGrowthRate

  double precision function adiabaticGnedin2004KSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the adiabaticGnedin2004 density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$), using the expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout)         :: self
    type            (treeNode                            ), intent(inout), target :: node
    double precision                                      , intent(in   )         :: waveNumber

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004KSpace=self%darkMatterProfileDMO_%kSpace         (node,waveNumber)
    else
       adiabaticGnedin2004KSpace=self                      %kSpaceNumerical(node,waveNumber)
    end if
    return
  end function adiabaticGnedin2004KSpace

  double precision function adiabaticGnedin2004FreefallRadius(self,node,time)
    !% Returns the freefall radius in the adiabaticGnedin2004 density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: time

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004FreefallRadius=self%darkMatterProfileDMO_%freefallRadius         (node,time)
    else
       adiabaticGnedin2004FreefallRadius=self                      %freefallRadiusNumerical(node,time)
    end if
    return
  end function adiabaticGnedin2004FreefallRadius

  double precision function adiabaticGnedin2004FreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the adiabaticGnedin2004 density profile at the specified
    !% {\normalfont \ttfamily time} (given in Gyr).
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: time

    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       adiabaticGnedin2004FreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate         (node,time)
    else
       adiabaticGnedin2004FreefallRadiusIncreaseRate=self                      %freefallRadiusIncreaseRateNumerical(node,time)
    end if
    return
  end function adiabaticGnedin2004FreefallRadiusIncreaseRate

  double precision function adiabaticGnedin2004RadiusInitial(self,node,radius)
    !% Compute the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    !% \cite{gnedin_response_2004}.
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    double precision                                      , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-2
    type            (rootFinder                          ), save          :: finder
    !$omp threadprivate(finder)
    integer                                                               :: i                      , j                       , &
         &                                                                   iMod
    double precision                                                      :: radiusUpperBound

    ! Reset stored solutions if the node has changed.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check for a previously computed solution.
    if (self%radiusPreviousIndexMaximum > 0 .and. any(self%radiusPrevious(1:self%radiusPreviousIndexMaximum) == radius)) then
       adiabaticGnedin2004RadiusInitial=0.0d0
       do i=1,self%radiusPreviousIndexMaximum
          if (self%radiusPrevious(i) == radius) then
             adiabaticGnedin2004RadiusInitial=self%radiusInitialPrevious(i)
             exit
          end if
       end do
       return
    end if
    ! Get the virial radius of the node.
    self%radiusVirial=self%darkMatterHaloScale_%virialRadius(node)
    ! Return radius unchanged if larger than the virial radius.
    if (radius >= self%radiusVirial) then
       adiabaticGnedin2004RadiusInitial=radius
       return
    end if
    ! Compute the various factors needed by this calculation.
    call self%computeFactors(node,radius,computeGradientFactors=.false.)
    ! If no baryons present at this radius, so initial radius is unchanged.
    if (self%baryonicFinalTerm <= 0.0d0) then
       adiabaticGnedin2004RadiusInitial=radius
       return
    end if
    ! Check that solution is within bounds.
    if (adiabaticGnedin2004Solver(self%radiusVirial) < 0.0d0) then
       adiabaticGnedin2004RadiusInitial=self%radiusVirial
       return
    end if
    j=-1
    if (self%radiusPreviousIndexMaximum > 0) then
       ! No exact match exists, look for approximate matches.
       do i=1,self%radiusPreviousIndexMaximum
          iMod=modulo(self%radiusPreviousIndex-i,adiabaticGnedin2004StoreCount)+1
          if (abs(radius-self%radiusPrevious(iMod))/self%radiusPrevious(iMod) < toleranceRelative) then
             j=iMod
             exit
          end if
       end do
    end if
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(adiabaticGnedin2004Solver                   )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
    end if
    ! Find the solution for initial radius.
    if (j == -1) then
       ! No previous solution to use as an initial guess. Instead, we make an estimate of the initial radius under the
       ! assumption that the mass of dark matter (in the initial profile) enclosed within the mean initial radius is the
       ! same as enclosed within the mean final radius. Since the initial and final radii are typically not too
       ! different, and since the mean radius is a weak (ѡ<1) function of the radius this is a useful
       ! approximation. Furthermore, since it will underestimate the actual mass within the initial mean radius it gives
       ! an overestimate of the initial radius. This means that we have a bracketing of the initial radius which we can
       ! use in the solver.
       radiusUpperBound   =  +(                                                                                        &
            &                  +self%baryonicFinalTerm                                                                 &
            &                  /self%darkMatterProfileDMO_%enclosedMass(node,self%radiusOrbitalMean(self%radiusFinal)) &
            &                  +self%darkMatterDistributedFraction                                                     &
            &                  *self%radiusFinal                                                                       &
            &                 )                                                                                        &
            &                /  self%initialMassFraction
       if (radiusUpperBound < radius) radiusUpperBound=radius
       call finder%rangeExpand(                                                             &
            &                  rangeExpandUpward            =1.1d0                        , &
            &                  rangeExpandDownward          =0.9d0                        , &
            &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                  rangeExpandType              =rangeExpandMultiplicative      &
            &                 )
       adiabaticGnedin2004RadiusInitial=finder%find(rootRange=[radius,radiusUpperBound])
    else
       ! Use previous solution as an initial guess.
       call finder%rangeExpand(                                                                   &
            &                  rangeExpandDownward          =1.0d0/sqrt(1.0d0+toleranceRelative), &
            &                  rangeExpandUpward            =1.0d0*sqrt(1.0d0+toleranceRelative), &
            &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative      , &
            &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive      , &
            &                  rangeExpandType              =rangeExpandMultiplicative            &
            &                 )
       adiabaticGnedin2004RadiusInitial=finder%find(                                                                        &
            &                                       rootRange=[                                                             &
            &                                                  self%radiusInitialPrevious(j)/sqrt(1.0d0+toleranceRelative), &
            &                                                  self%radiusInitialPrevious(j)*sqrt(1.0d0+toleranceRelative)  &
            &                                                 ]                                                             &
            &                                      )
    end if
    ! Store this solution.
    self%radiusPreviousIndex                                 =modulo(self%radiusPreviousIndex         ,adiabaticGnedin2004StoreCount)+1
    self%radiusPreviousIndexMaximum                          =min   (self%radiusPreviousIndexMaximum+1,adiabaticGnedin2004StoreCount)
    self%radiusPrevious            (self%radiusPreviousIndex)=radius
    self%radiusInitialPrevious     (self%radiusPreviousIndex)=adiabaticGnedin2004RadiusInitial
    return
  end function adiabaticGnedin2004RadiusInitial

  double precision function adiabaticGnedin2004RadiusInitialDerivative(self,node,radius)
    !% Compute the derivative of the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    !% \cite{gnedin_response_2004}.
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: radius
    double precision                                      , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder                          ), save          :: finder
    !$omp threadprivate(finder)

    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rangeExpand (                                                             &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeDownwardLimit           =0.0d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
       call finder%rootFunction(adiabaticGnedin2004DerivativeSolver                  )
       call finder%tolerance   (toleranceAbsolute                  ,toleranceRelative)
    end if
    ! Reset stored solutions if the node has changed.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Compute the various factors needed by this calculation.
    call self%computeFactors(node,radius,computeGradientFactors=.true.)
    ! Return unit derivative if radius is larger than the virial radius.
    if (radius >= self%radiusVirial) then
       adiabaticGnedin2004RadiusInitialDerivative=1.0d0
       return
    end if
    ! Compute initial radius, and derivatives of initial and final mean radii.
    self%radiusFinal                    =                                           radius
    self%radiusInitial_                 =self%radiusInitial              (node,     radius        )
    self%radiusInitialMeanSelfDerivative=self%radiusOrbitalMeanDerivative(     self%radiusInitial_)
    self%radiusFinalMeanSelfDerivative  =self%radiusOrbitalMeanDerivative(          radius        )
    ! Find the solution for initial radius.
    adiabaticGnedin2004RadiusInitialDerivative=finder%find(rootGuess=1.0d0)
    return
  end function adiabaticGnedin2004RadiusInitialDerivative

  subroutine adiabaticGnedin2004ComputeFactors(self,node,radius,computeGradientFactors)
    !% Compute various factors needed when solving for the initial radius in the dark matter halo using the adiabatic contraction
    !% algorithm of \cite{gnedin_response_2004}.
    use :: Galacticus_Nodes            , only : nodeComponentBasic             , optimizeForEnclosedMassSummation, optimizeForRotationCurveGradientSummation, optimizeForRotationCurveSummation, &
          &                                     reductionSummation             , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    !# <include directive="rotationCurveTask" name="radiusSolverRotationCurveTask" type="moduleUse">
    !# <exclude>Dark_Matter_Profile_Structure_Tasks</exclude>
    include 'dark_matter_profiles.nonDMO.adiabatic_Gnedin2004.rotation_curve.tasks.modules.inc'
    !# </include>
    !# <include directive="rotationCurveGradientTask" name="radiusSolverRotationCurveGradientTask" type="moduleUse">
    !# <exclude>Dark_Matter_Profile_Structure_Tasks</exclude>
    include 'dark_matter_profiles.nonDMO.adiabatic_Gnedin2004.rotation_curve_gradient.tasks.modules.inc'
    !# </include>
    !# <include directive="enclosedMassTask" name="radiusSolverEnclosedMassTask" type="moduleUse">
    !# <exclude>Dark_Matter_Profile_Structure_Tasks</exclude>
    include 'dark_matter_profiles.nonDMO.adiabatic_Gnedin2004.enclosed_mass.tasks.modules.inc'
    !# </include>
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout), target  :: self
    type            (treeNode                            ), intent(inout), target  :: node
    double precision                                      , intent(in   )          :: radius
    logical                                               , intent(in   )          :: computeGradientFactors
    type            (treeNode                            )               , pointer :: nodeCurrent                           , nodeHost
    class           (nodeComponentBasic                  )               , pointer :: basic
    double precision                                      , parameter              :: toleranceAbsolute               =0.0d0, toleranceRelative   =1.0d-3
    procedure       (adiabaticGnedin2004MassEnclosed     )               , pointer :: mapFunction
    double precision                                                               :: massBaryonicSelfTotal                 , massBaryonicTotal          , &
         &                                                                            massComponent                         , velocityComponent          , &
         &                                                                            velocityComponentSquaredGradient      , velocityCircularSquared    , &
         &                                                                            velocityCircularSquaredGradient

    ! Set module-scope pointers to node and self.
    adiabaticGnedin2004Node => node
    adiabaticGnedin2004Self => self
    ! Store the final radius and its orbit-averaged mean.
    self%radiusFinal    =                       radius
    self%radiusFinalMean=self%radiusOrbitalMean(radius)
    self%radiusShared   =self%radiusFinalMean
    ! Compute the baryonic contribution to the rotation curve.
    mapFunction             => adiabaticGnedin2004VelocityCircularSquared
    velocityCircularSquared =  node%mapDouble0(mapFunction,reductionSummation,optimizeFor=optimizeForRotationCurveSummation)
    !# <include directive="rotationCurveTask" name="radiusSolverRotationCurveTask" type="functionCall" functionType="function" returnParameter="velocityComponent">
    !#  <exclude>Dark_Matter_Profile_Rotation_Curve_Task</exclude>
    !#  <functionArgs>node,self%radiusFinalMean,adiabaticGnedin2004ComponentType,adiabaticGnedin2004MassType</functionArgs>
    !#  <onReturn>velocityCircularSquared=velocityCircularSquared+velocityComponent**2</onReturn>
    include 'dark_matter_profiles.nonDMO.adiabatic_Gnedin2004.rotation_curve.tasks.inc'
    !# </include>
    self%baryonicFinalTerm=velocityCircularSquared*self%radiusFinalMean*self%radiusFinal/gravitationalConstantGalacticus
    ! Compute the baryonic contribution to the rotation curve.
    if (computeGradientFactors) then
       mapFunction                     => adiabaticGnedin2004VelocityCircularSquaredGradient
       velocityCircularSquaredGradient =  node%mapDouble0(mapFunction,reductionSummation,optimizeFor=optimizeForRotationCurveGradientSummation)
       !# <include directive="rotationCurveGradientTask" name="radiusSolverRotationCurveGradientTask" type="functionCall" functionType="function" returnParameter="velocityComponentSquaredGradient">
       !#  <exclude>Dark_Matter_Profile_Rotation_Curve_Gradient_Task</exclude>
       !#  <functionArgs>node,self%radiusFinalMean,adiabaticGnedin2004ComponentType,adiabaticGnedin2004MassType</functionArgs>
       !#  <onReturn>velocityCircularSquaredGradient=velocityCircularSquaredGradient+velocityComponentSquaredGradient</onReturn>
       include 'dark_matter_profiles.nonDMO.adiabatic_Gnedin2004.rotation_curve_gradient.tasks.inc'
       !# </include>
       self%baryonicFinalTermDerivative=velocityCircularSquaredGradient*self%radiusOrbitalMeanDerivative(self%radiusFinal)*self%radiusFinalMean*self%radiusFinal/gravitationalConstantGalacticus
    end if
    ! Compute the initial baryonic contribution from this halo, and any satellites.
    if (.not.self%massesComputed) then
       massBaryonicTotal     =  0.0d0
       massBaryonicSelfTotal =  0.0d0
       nodeCurrent           => node
       nodeHost              => node
       do while (associated(nodeCurrent))
          mapFunction       => adiabaticGnedin2004MassEnclosed
          massBaryonicTotal =  massBaryonicTotal+nodeCurrent%mapDouble0(mapFunction,reductionSummation,optimizeFor=optimizeForEnclosedMassSummation)
          !# <include directive="enclosedMassTask" name="radiusSolverEnclosedMassTask" type="functionCall" functionType="function" returnParameter="massComponent">
          !#  <exclude>Dark_Matter_Profile_Enclosed_Mass_Task</exclude>
          !#  <functionArgs>nodeCurrent,radiusLarge,adiabaticGnedin2004ComponentType,adiabaticGnedin2004MassType,adiabaticGnedin2004WeightBy,adiabaticGnedin2004WeightIndex</functionArgs>
          !#  <onReturn>massBaryonicTotal=massBaryonicTotal+massComponent</onReturn>
          include 'dark_matter_profiles.nonDMO.adiabatic_Gnedin2004.enclosed_mass.tasks.inc'
          !# </include>
          if (associated(nodeCurrent,nodeHost)) then
             massBaryonicSelfTotal=massBaryonicTotal
             do while (associated(nodeCurrent%firstSatellite))
                nodeCurrent => nodeCurrent%firstSatellite
             end do
             if (associated(nodeCurrent,nodeHost)) nodeCurrent => null()
          else
             if (associated(nodeCurrent%sibling)) then
                nodeCurrent => nodeCurrent%sibling
                do while (associated(nodeCurrent%firstSatellite))
                   nodeCurrent => nodeCurrent%firstSatellite
                end do
             else
                nodeCurrent => nodeCurrent%parent
                if (associated(nodeCurrent,node)) nodeCurrent => null()
             end if
          end if
       end do
       ! Limit masses to physical values.
       massBaryonicSelfTotal=max(massBaryonicSelfTotal,0.0d0)
       massBaryonicTotal    =max(massBaryonicTotal    ,0.0d0)
       ! Compute the fraction of matter assumed to be distributed like the dark matter.
       basic                              => node%basic()
       self%darkMatterDistributedFraction =min((self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon())/self%cosmologyParameters_%OmegaMatter()+(massBaryonicTotal-massBaryonicSelfTotal)/basic%mass(),1.0d0)
       ! Compute the initial mass fraction.
       self%initialMassFraction           =min((self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon())/self%cosmologyParameters_%OmegaMatter()+ massBaryonicTotal                       /basic%mass(),1.0d0)
       ! Record that masses (and mass fractions) have been computed.
       self%massesComputed=.true.
    end if
    return
  end subroutine adiabaticGnedin2004ComputeFactors

  double precision function adiabaticGnedin2004RadiusOrbitalMean(self,radius)
    !% Returns the orbit averaged radius for dark matter corresponding the given {\normalfont \ttfamily radius} using the model of
    !% \cite{gnedin_response_2004}.
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius

    adiabaticGnedin2004RadiusOrbitalMean=+self%A                                                          &
         &                               *self%radiusVirial                                               &
         &                               *self%radiusExponentiator%exponentiate(radius/self%radiusVirial)
    return
  end function adiabaticGnedin2004RadiusOrbitalMean

  double precision function adiabaticGnedin2004RadiusOrbitalMeanDerivative(self,radius)
    !% Returns the derivative of the orbit averaged radius for dark matter corresponding the given {\normalfont \ttfamily radius} using the model of
    !% \cite{gnedin_response_2004}.
    implicit none
    class           (darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius

    adiabaticGnedin2004RadiusOrbitalMeanDerivative=+self%A                &
         &                                         *(                     &
         &                                           +     radius         &
         &                                           /self%radiusVirial   &
         &                                          )**(self%omega-1.0d0)
    return
  end function adiabaticGnedin2004RadiusOrbitalMeanDerivative

  double precision function adiabaticGnedin2004MassEnclosed(component)
    !% Unary function returning the enclosed mass in a component. Suitable for mapping over components. Ignores the dark matter
    !% profile.
    use :: Galacticus_Nodes, only : nodeComponent, nodeComponentDarkMatterProfile
    implicit none
    class(nodeComponent), intent(inout) :: component

    select type (component)
    class is (nodeComponentDarkMatterProfile)
       adiabaticGnedin2004MassEnclosed=0.0d0
    class default
       adiabaticGnedin2004MassEnclosed=component%enclosedMass(radiusLarge,adiabaticGnedin2004ComponentType,adiabaticGnedin2004MassType,adiabaticGnedin2004WeightBy,adiabaticGnedin2004WeightIndex)
    end select
    return
  end function adiabaticGnedin2004MassEnclosed

  double precision function adiabaticGnedin2004VelocityCircularSquared(component)
    !% Unary function returning the squared rotation curve in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent, nodeComponentDarkMatterProfile
    implicit none
    class(nodeComponent), intent(inout) :: component

    select type (component)
    class is (nodeComponentDarkMatterProfile)
       adiabaticGnedin2004VelocityCircularSquared=0.0d0
    class default
       adiabaticGnedin2004VelocityCircularSquared=component%rotationCurve(adiabaticGnedin2004Self%radiusShared,adiabaticGnedin2004ComponentType,adiabaticGnedin2004MassType)**2
    end select
    return
  end function adiabaticGnedin2004VelocityCircularSquared

  double precision function adiabaticGnedin2004VelocityCircularSquaredGradient(component)
    !% Unary function returning the squared rotation curve gradient in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent, nodeComponentDarkMatterProfile
    implicit none
    class(nodeComponent), intent(inout) :: component

    select type (component)
    class is (nodeComponentDarkMatterProfile)
       adiabaticGnedin2004VelocityCircularSquaredGradient=0.0d0
    class default
       adiabaticGnedin2004VelocityCircularSquaredGradient=component%rotationCurveGradient(adiabaticGnedin2004Self%radiusShared,adiabaticGnedin2004ComponentType,adiabaticGnedin2004MassType)
    end select
    return
  end function adiabaticGnedin2004VelocityCircularSquaredGradient

  double precision function adiabaticGnedin2004Solver(radiusInitial)
    !% Root function used in finding the initial radius in the dark matter halo when solving for adiabatic contraction.
    implicit none
    double precision, intent(in   ) :: radiusInitial
    double precision                :: massDarkMatterInitial, radiusInitialMean

    ! Find the initial mean orbital radius.
    radiusInitialMean    =adiabaticGnedin2004Self                      %radiusOrbitalMean(                        radiusInitial    )
    ! Get the mass of dark matter inside the initial radius.
    massDarkMatterInitial=adiabaticGnedin2004Self%darkMatterProfileDMO_%enclosedMass     (adiabaticGnedin2004Node,radiusInitialMean)
    ! Compute the root function.
    adiabaticGnedin2004Solver=+massDarkMatterInitial                                                                          &
         &                    *(                                                                                              &
         &                      +adiabaticGnedin2004Self%initialMassFraction*                                   radiusInitial &
         &                      -adiabaticGnedin2004Self%darkMatterDistributedFraction *adiabaticGnedin2004Self%radiusFinal   &
         &                     )                                                                                              &
         &                    -adiabaticGnedin2004Self%baryonicFinalTerm
    return
  end function adiabaticGnedin2004Solver

  double precision function adiabaticGnedin2004DerivativeSolver(radiusInitialDerivative)
    !% Root function used in finding the derivative of the initial radius in the dark matter halo when solving for adiabatic
    !% contraction.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radiusInitialDerivative
    double precision                :: densityDarkMatterInitial, massDarkMatterInitial, &
         &                             radiusInitialMean

    ! Find the initial mean orbital radius.
    radiusInitialMean       =adiabaticGnedin2004Self                      %radiusOrbitalMean(                        adiabaticGnedin2004Self%radiusInitial_   )
    ! Get the mass of dark matter inside the initial radius.
    massDarkMatterInitial   =adiabaticGnedin2004Self%darkMatterProfileDMO_%enclosedMass     (adiabaticGnedin2004Node,                        radiusInitialMean)
    ! Get the mass of dark matter inside the initial radius.
    densityDarkMatterInitial=adiabaticGnedin2004Self%darkMatterProfileDMO_%density          (adiabaticGnedin2004Node,                        radiusInitialMean)
    ! Compute the root function.
    adiabaticGnedin2004DerivativeSolver=+massDarkMatterInitial                                                                                   &
         &                              *(                                                                                                       &
         &                                +adiabaticGnedin2004Self%initialMassFraction*                                  radiusInitialDerivative &
         &                                -adiabaticGnedin2004Self%darkMatterDistributedFraction                                                 &
         &                               )                                                                                                       &
         &                              +(                                                                                                       &
         &                                +adiabaticGnedin2004Self%initialMassFraction          *adiabaticGnedin2004Self%radiusInitial_          &
         &                                -adiabaticGnedin2004Self%darkMatterDistributedFraction*adiabaticGnedin2004Self%radiusFinal             &
         &                               )                                                                                                       &
         &                              *4.0d0                                                                                                   &
         &                              *Pi                                                                                                      &
         &                              *radiusInitialMean**2                                                                                    &
         &                              *densityDarkMatterInitial                                                                                &
         &                              *adiabaticGnedin2004Self%radiusInitialMeanSelfDerivative                                                 &
         &                              *               radiusInitialDerivative                                                                  &
         &                              -adiabaticGnedin2004Self%baryonicFinalTerm                                                               &
         &                              *(                                                                                                       &
         &                                +1.0d0                                                /adiabaticGnedin2004Self%radiusFinal             &
         &                                +adiabaticGnedin2004Self%radiusFinalMeanSelfDerivative/adiabaticGnedin2004Self%radiusFinalMean         &
         &                               )                                                                                                       &
         &                              -adiabaticGnedin2004Self%baryonicFinalTermDerivative
    return
  end function adiabaticGnedin2004DerivativeSolver
