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

  !!{RST
  Implements the :cite:t:`gnedin_effect_2000` filtering mass calculation.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctions      , cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParameters     , cosmologyParametersClass
  use :: Intergalactic_Medium_State, only : intergalacticMediumState, intergalacticMediumStateClass
  use :: Linear_Growth             , only : linearGrowth            , linearGrowthClass
  use :: Tables                    , only : table                   , table1D

  public :: gnedin2000ODEs

  !![
  <intergalacticMediumFilteringMass name="intergalacticMediumFilteringMassGnedin2000" docformat="rst">
   <description>
   An implementation of the :cite:t:`gnedin_effect_2000` filtering mass calculation, which determines the characteristic halo mass below which gas accretion is suppressed by photoionization heating from the intergalactic radiation field. The filtering mass is computed by integrating an ODE system driven by the :term:`IGM` thermal state and linear growth history.
   </description>
  </intergalacticMediumFilteringMass>
  !!]
  type, extends(intergalacticMediumFilteringMassClass) :: intergalacticMediumFilteringMassGnedin2000
     !!{RST
     An implementation of the :cite:t:`gnedin_effect_2000` filtering mass calculation.
     !!}
     private
     class           (cosmologyParametersClass     ), pointer     :: cosmologyParameters_      => null()
     class           (cosmologyFunctionsClass      ), pointer     :: cosmologyFunctions_       => null()
     class           (intergalacticMediumStateClass), pointer     :: intergalacticMediumState_ => null()
     class           (linearGrowthClass            ), pointer     :: linearGrowth_             => null()
     logical                                                      :: initialized
     integer                                                      :: countTimes
     double precision                                             :: timeMaximum                        , timeMinimum
     logical                                                      :: timeTooEarlyIsFatal
     class           (table1D                      ), allocatable :: table
     type            (varying_string               )              :: fileName
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulate the filtering mass to encompass at least the given ``time``." method="tabulate" />
       <method description="Set the initial conditions for the ODE system." method="conditionsInitialODEs" />
       <method description="Return coefficients for the early-epoch fitting function to the filtering mass." method="coefficientsEarlyEpoch" />
       <method description="Return the early-epoch solution for the filtering mass." method="massFilteringEarlyEpoch" />
     </methods>
     !!]
     final     ::                                gnedin2000Destructor
     procedure :: massFiltering               => gnedin2000MassFiltering
     procedure :: massFilteringRateOfChange   => gnedin2000MassFilteringRateOfChange
     procedure :: fractionBaryons             => gnedin2000FractionBaryons
     procedure :: fractionBaryonsRateOfChange => gnedin2000FractionBaryonsRateOfChange
     procedure :: fractionBaryonsGradientMass => gnedin2000FractionBaryonsGradientMass
     procedure :: tabulate                    => gnedin2000Tabulate
     procedure :: conditionsInitialODEs       => gnedin2000ConditionsInitialODEs
     procedure :: coefficientsEarlyEpoch      => gnedin2000CoefficientsEarlyEpoch
     procedure :: massFilteringEarlyEpoch     => gnedin2000MassFilteringEarlyEpoch
  end type intergalacticMediumFilteringMassGnedin2000

  interface intergalacticMediumFilteringMassGnedin2000
     !!{RST
     Constructors for the filtering mass class.
     !!}
     module procedure gnedin2000ConstructorParameters
     module procedure gnedin2000ConstructorInternal
  end interface intergalacticMediumFilteringMassGnedin2000

  ! Parameter controlling fine-grainedness of filtering mass tabulations.
  integer                , parameter :: tablePointsPerDecade=100

  ! Interval, measured in lattice steps, to which the bounds of the tabulated range are pinned. Half decades are used rather than
  ! whole decades since a whole decade of cosmic time spans most of the history of the universe.
  integer                , parameter :: tableAnchorEvery    =tablePointsPerDecade/2

contains

  function gnedin2000ConstructorParameters(parameters) result(self)
    !!{RST
    Default constructor for the file :term:`IGM` state class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (intergalacticMediumFilteringMassGnedin2000)                :: self
    type   (inputParameters                           ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_
    class  (cosmologyParametersClass                  ), pointer       :: cosmologyParameters_
    class  (linearGrowthClass                         ), pointer       :: linearGrowth_
    class  (intergalacticMediumStateClass             ), pointer       :: intergalacticMediumState_
    logical                                                            :: timeTooEarlyIsFatal

    !![
    <inputParameter docformat="rst">
      <name>timeTooEarlyIsFatal</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, requesting the filtering mass at a time earlier than the initial time provided by the :cite:t:`naoz_growth_2005` fit will result in a fatal error. Otherwise, the filtering mass is fixed at this initial value for earlier times.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="intergalacticMediumState" name="intergalacticMediumState_" source="parameters"/>
    !!]
    self=intergalacticMediumFilteringMassGnedin2000(timeTooEarlyIsFatal,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,intergalacticMediumState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="intergalacticMediumState_"/>
    !!]
    return
  end function gnedin2000ConstructorParameters

  function gnedin2000ConstructorInternal(timeTooEarlyIsFatal,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,intergalacticMediumState_) result(self)
    !!{RST
    Constructor for the filtering mass class.
    !!}
    use :: File_Utilities    , only : Directory_Make      , File_Path
    use :: ISO_Varying_String, only : char
    use :: Table_Caches      , only : Table_Cache_File_Name
    implicit none
    type   (intergalacticMediumFilteringMassGnedin2000)                        :: self
    class  (cosmologyParametersClass                  ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass                   ), intent(in   ), target :: cosmologyFunctions_
    class  (linearGrowthClass                         ), intent(in   ), target :: linearGrowth_
    class  (intergalacticMediumStateClass             ), intent(in   ), target :: intergalacticMediumState_
    logical                                            , intent(in   )         :: timeTooEarlyIsFatal
    !![
    <constructorAssign variables="timeTooEarlyIsFatal, *cosmologyParameters_, *cosmologyFunctions_, *linearGrowth_, *intergalacticMediumState_"/>
    !!]

    self%initialized=.false.
    self%fileName   =Table_Cache_File_Name(                                                                                   &
         &                                 subDirectory    ='intergalacticMedium'                                           , &
         &                                 objectType      =char(self%objectType      (                                   )), &
         &                                 hashedDescriptor=char(self%hashedDescriptor(includeSourceDigest         =.true.  , &
         &                                                                             includeFileModificationTimes=.true.))  &
         &                                )
    call Directory_Make(File_Path(self%fileName))
    return
  end function gnedin2000ConstructorInternal

  subroutine gnedin2000Destructor(self)
    !!{RST
    Destructor for the filtering mass class.
    !!}
    implicit none
    type (intergalacticMediumFilteringMassGnedin2000), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%intergalacticMediumState_"/>
    !!]
    return
  end subroutine gnedin2000Destructor

  double precision function gnedin2000MassFiltering(self,time)
    !!{RST
    Return the filtering mass at the given ``time``.
    !!}
    implicit none
    class           (intergalacticMediumFilteringMassGnedin2000), intent(inout) :: self
    double precision                                            , intent(in   ) :: time

    call self%tabulate(time)
    gnedin2000MassFiltering=self%table%interpolate(max(time,self%timeMinimum))
    return
  end function gnedin2000MassFiltering

  double precision function gnedin2000MassFilteringRateOfChange(self,time)
    !!{RST
    Return the rate of change of the filtering mass at the given ``time``.
    !!}
    implicit none
    class           (intergalacticMediumFilteringMassGnedin2000), intent(inout) :: self
    double precision                                            , intent(in   ) :: time

    call self%tabulate(time)
    gnedin2000MassFilteringRateOfChange=self%table%interpolateGradient(max(time,self%timeMinimum))
    return
  end function gnedin2000MassFilteringRateOfChange

  double precision function gnedin2000FractionBaryons(self,mass,time)
    !!{RST
    Return the rate of change of the fraction of baryons accreted into a halo of the given ``mass`` at the ``time``.
    !!}
    implicit none
    class           (intergalacticMediumFilteringMassGnedin2000), intent(inout) :: self
    double precision                                            , intent(in   ) :: mass, time

    gnedin2000FractionBaryons=1.0d0/(1.0d0+(2.0d0**(1.0d0/3.0d0)-1.0d0)*8.0d0*self%massFiltering(time)/mass)**3
    return
  end function gnedin2000FractionBaryons

  double precision function gnedin2000FractionBaryonsRateOfChange(self,mass,time)
    !!{RST
    Return the rate of change of the fraction of baryons accreted into a halo of the given ``mass`` at the ``time``.
    !!}
    implicit none
    class           (intergalacticMediumFilteringMassGnedin2000), intent(inout) :: self
    double precision                                            , intent(in   ) :: mass, time

    gnedin2000FractionBaryonsRateOfChange=-3.0d0                                                    &
         &                                *(2.0d0**(1.0d0/3.0d0)-1.0d0)*8.0d0                       &
         &                                *self%massFilteringRateOfChange(     time)                &
         &                                /     mass                                                &
         &                                *self%fractionBaryons          (mass,time)**(4.0d0/3.0d0)
    return
  end function gnedin2000FractionBaryonsRateOfChange

  double precision function gnedin2000FractionBaryonsGradientMass(self,mass,time)
    !!{RST
    Return the gradient with respect to mass of the fraction of baryons accreted into a halo of the given ``mass`` at the ``time``.
    !!}
    implicit none
    class           (intergalacticMediumFilteringMassGnedin2000), intent(inout) :: self
    double precision                                            , intent(in   ) :: mass, time

    gnedin2000FractionBaryonsGradientMass=-(2.0d0**(1.0d0/3.0d0)-1.0d0)*8.0d0             &
         &                                *self%massFiltering  (     time)                &
         &                                /     mass                      **2             &
         &                                *self%fractionBaryons(mass,time)**(4.0d0/3.0d0)
    return
  end function gnedin2000FractionBaryonsGradientMass

  subroutine gnedin2000Tabulate(self,time)
    !!{RST
    Construct a table of filtering mass as a function of cosmological time.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_ODE_Solvers   , only : odeSolver
    use :: Numerical_Ranges        , only : Range_Pinned            , rangeLattice     , gridSchemePerDecade
    use :: Table_Caches            , only : Table_Cache_Restore     , Table_Cache_Store
    use :: Tables                  , only : table1DLogarithmicLinear
    implicit none
    class           (intergalacticMediumFilteringMassGnedin2000), intent(inout), target :: self
    double precision                                            , intent(in   )         :: time
    double precision                                            , parameter             :: redshiftMaximumNaozBarkana=150.0d+0 ! Maximum redshift at which fitting function of Naoz & Barkana is valid.
    double precision                                            , dimension(3)          :: massFiltering                      , massFilteringScales
    double precision                                            , parameter             :: odeToleranceAbsolute      =  1.0d-3, odeToleranceRelative=1.0d-3
    type            (odeSolver                                 )                        :: solver
    type            (rangeLattice                              )                        :: lattice
    logical                                     , allocatable   , dimension(:)          :: isComputed
    integer                                                                             :: iTime                              , status
    double precision                                                                    :: timeInitial                        , timeCurrent                , &
         &                                                                                 timePresent

    ! Evaluate a suitable starting time for filtering mass calculations. The fitting function of Naoz & Barkana used to set the
    ! initial conditions is valid only at redshifts below 150, which therefore imposes a hard lower limit on the range which can
    ! be tabulated.
    timeInitial=self%cosmologyFunctions_ %cosmicTime                 (                            &
         &       self%cosmologyFunctions_%expansionFactorFromRedshift (                           &
         &                                                             redshiftMaximumNaozBarkana &
         &                                                            )                           &
         &                                                           )
    timePresent=self%cosmologyFunctions_ %cosmicTime                 (1.0d0)
    ! Abort if time is too early.
    if (time <= timeInitial .and. self%timeTooEarlyIsFatal) call Error_Report('time is too early'//{introspection:location})
    if (.not.allocated(self%table)) allocate(table1DLogarithmicLinear :: self%table)
    ! Find the range of times to tabulate, pinned to an absolute lattice of points per decade so that the tabulation is
    ! independent of the time at which it was first requested. The table always spans the present epoch, and is never tabulated
    ! beyond it unless a later time is requested.
    lattice=Range_Pinned(                                          &
         &                              time                     , &
         &                              tablePointsPerDecade     , &
         &                              gridSchemePerDecade      , &
         &               rangeCurrent  =[timePresent,timePresent], &
         &               latticeCurrent=self%table%lattice       , &
         &               limitMinimum  =timeInitial              , &
         &               limitMaximum  =max(timePresent,time)    , &
         &               anchorEvery   =tableAnchorEvery           &
         &              )
    if (self%table%lattice%covers(lattice)) return
    ! Merge in any tabulation already cached on disk, then re-evaluate the range required in the light of it.
    call Table_Cache_Restore(self%table,self%fileName,status)
    lattice=Range_Pinned(                                          &
         &                              time                     , &
         &                              tablePointsPerDecade     , &
         &                              gridSchemePerDecade      , &
         &               rangeCurrent  =[timePresent,timePresent], &
         &               latticeCurrent=self%table%lattice       , &
         &               limitMinimum  =timeInitial              , &
         &               limitMaximum  =max(timePresent,time)    , &
         &               anchorEvery   =tableAnchorEvery           &
         &              )
    call self%table%extend(lattice,isComputed)
    ! Solve the ODE system for those points which are not already known. Note that no lock is held while doing so.
    if (any(.not.isComputed)) then
       ! Set the initial state for the composite variables used to solve for filtering mass - this also establishes the scales
       ! used by the ODE solver, so must be done before it is constructed.
       call self%conditionsInitialODEs(timeInitial,massFiltering,massFilteringScales)
       solver=odeSolver(3_c_size_t,massFilteringODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative,scale=massFilteringScales)
       select type (table_ => self%table)
       type is (table1DLogarithmicLinear)
          do iTime=1,lattice%count
             if (isComputed(iTime)) cycle
             ! Reset the initial state composite variables used to solve for filtering mass.
             call self%conditionsInitialODEs(timeInitial,massFiltering,massFilteringScales)
             ! Solve the ODE system
             timeCurrent=timeInitial
             if (table_%x(iTime) > timeInitial) &
                  & call solver%solve(timeCurrent,table_%x(iTime),massFiltering)
             call table_%populate(massFiltering(3),iTime)
          end do
       end select
       ! Merge with, and write back, the cache file.
       call Table_Cache_Store(self%table,self%fileName)
    end if
    ! Record the tabulated range.
    self%timeMinimum=self%table%lattice%minimum()
    self%timeMaximum=self%table%lattice%maximum()
    self%countTimes =self%table%lattice%count
    self%initialized=.true.
    return

  contains

    integer function massFilteringODEs(time,properties,propertiesRateOfChange)
      !!{RST
      Evaluates the ODEs controlling the evolution temperature.
      !!}
      use :: Interface_GSL                   , only : GSL_Success
      use :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial, hydrogenByMassPrimordial
      use :: Numerical_Constants_Atomic      , only : electronMass          , massHeliumAtom          , massHydrogenAtom
      implicit none
      double precision, intent(in  )                :: time
      double precision, intent(in   ), dimension(:) :: properties
      double precision, intent(  out), dimension(:) :: propertiesRateOfChange
      double precision                              :: temperature            , massParticleMean

       ! Find mean particle mass.
      massParticleMean=+(hydrogenByMassPrimordial*(1.0d0+self%intergalacticMediumState_%electronFraction(time)*electronMass/massHydrogenAtom)                 +heliumByMassPrimordial               ) &
           &           /(hydrogenByMassPrimordial*(1.0d0+self%intergalacticMediumState_%electronFraction(time)                              )/massHydrogenAtom+heliumByMassPrimordial/massHeliumAtom)
      ! Get the temperature.
      temperature=self%intergalacticMediumState_%temperature(time)
      ! Compute the rates of change for the ODE system.
      propertiesRateOfChange=gnedin2000ODEs(self%cosmologyParameters_,self%cosmologyFunctions_,self%linearGrowth_,time,massParticleMean,temperature,properties)
      ! Return success.
      massFilteringODEs=GSL_Success
      return
    end function massFilteringODEs

  end subroutine gnedin2000Tabulate

  function gnedin2000MassFilteringEarlyEpoch(self,time) result (massFiltering)
    !!{RST
    Fitting function for the filtering mass at early epochs from :cite:t:`naoz_formation_2007`. Checks for valid range of redshift and cosmology for the fit to be valid.
    !!}
    implicit none
    class           (intergalacticMediumFilteringMassGnedin2000), intent(inout) :: self
    double precision                                                            :: massFiltering
    double precision                                            , intent(in   ) :: time
    double precision                                            , dimension(4)  :: coefficients
    double precision                                                            :: expansionFactor

    ! Compute the expansion factor.
    expansionFactor=self%cosmologyFunctions_%expansionFactor        (time)
    ! Get fitting function coefficients.
    coefficients   =self                    %coefficientsEarlyEpoch (time)
    ! Evaluate fitting function.
    massFiltering  =+exp(                                            &
         &               +coefficients(1)*(-log(expansionFactor))**3 &
         &               +coefficients(2)*(-log(expansionFactor))**2 &
         &               +coefficients(3)*(-log(expansionFactor))    &
         &               +coefficients(4)                            &
         &              )
    ! Compute the
    return
  end function gnedin2000MassFilteringEarlyEpoch

  subroutine gnedin2000ConditionsInitialODEs(self,time,massFilteringODEs,massFilteringScales)
    !!{RST
    Compute initial conditions for a system of three variables used to solve for the evolution of the filtering mass. The ODE system to be solved is

    .. math::

       \dot{y}_1 &=& y_2 \\
       \dot{y}_2 &=& -2 (\dot{a}/a) y_2 + (1+r_\mathrm{LSS}(t)) f_\mathrm{DM} D(t) \mathrm{k}_\mathrm{B} T(t)/\mu m_\mathrm{H} a^{-2} \\
       \dot{y}_3 &=& - 4 \pi^4 \bar{\rho} \dot{k}_\mathrm{F}(t)/ k_\mathrm{F}^4(t)

    with initial conditions

    .. math::

       y_1 &=& D(t)/k_\mathrm{F}^2(t) \\
       y_2 &=& \dot{y}_1 \\
       y_3 &=& M_\mathrm{F}(t)

    and where

    .. math::

       k_\mathrm{F}(t) = \pi / [M_\mathrm{F}(t) 3 / 4 \pi \bar{\rho}]^{1/3}

    , and :math:`r_\mathrm{LSS}(t)` is the function defined by :cite:t:`naoz_formation_2007`.
    !!}
    use :: Cosmology_Parameters    , only : hubbleUnitsTime
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (intergalacticMediumFilteringMassGnedin2000), intent(inout)                         :: self
    double precision                                            , intent(in   )                         :: time
    double precision                                            , intent(  out), dimension(3)           :: massFilteringODEs
    double precision                                            , intent(  out), dimension(3), optional :: massFilteringScales
    double precision                                                           , dimension(4)           :: coefficients
    double precision                                                                                    :: expansionFactor     , expansionRate      , &
         &                                                                                                 massFiltering       , wavenumberFiltering

    ! Get fitting function coefficients.
    coefficients          =self                    %coefficientsEarlyEpoch (time           )
    ! Compute expansion factor and rate.
    expansionFactor       =self%cosmologyFunctions_%expansionFactor        (time           )
    expansionRate         =self%cosmologyFunctions_%expansionRate          (expansionFactor)
    ! Get the filtering mass at the initial time.
    massFiltering         =self                    %massFilteringEarlyEpoch(time           )
    ! Find the corresponding filtering wavenumber.
    wavenumberFiltering   =+Pi                                            &
         &                 /(                                             &
         &                   +massFiltering                               &
         &                   *3.0d0                                       &
         &                   /4.0d0                                       &
         &                   /Pi                                          &
         &                   /self%cosmologyParameters_%OmegaMatter    () &
         &                   /self%cosmologyParameters_%densityCritical() &
         &                  )**(1.0d0/3.0d0)
    ! Evaluate the three ODE variables at the initial time.
    massFilteringODEs  (1)=+self%linearGrowth_%value                                 (time) &
         &                 /wavenumberFiltering**2
    massFilteringODEs  (2)=+(                                                               &
         &                   +self%linearGrowth_%value                               (time) &
         &                   *self%linearGrowth_%logarithmicDerivativeExpansionFactor(time) &
         &                   /wavenumberFiltering**2                                        &
         &                   +2.0d0                                                         &
         &                   /3.0d0                                                         &
         &                   *self%linearGrowth_%value                               (time) &
         &                   /wavenumberFiltering**2                                        &
         &                   *(                                                             &
         &                     -3.0d0                                                       &
         &                     *coefficients(1)                                             &
         &                     *log(expansionFactor)**2                                     &
         &                     +2.0d0                                                       &
         &                     *coefficients(2)                                             &
         &                     *log(expansionFactor)                                        &
         &                     -coefficients(3)                                             &
         &                    )                                                             &
         &                  )                                                               &
         &                 *expansionRate
    massFilteringODEs  (3)=+massFiltering
    ! Evaluate suitable absolute tolerance scales for the ODE variables.
    if (present(massFilteringScales)) then
       massFilteringScales(1)=+self%linearGrowth_%value(time)                               &
            &                 /wavenumberFiltering**2
       massFilteringScales(2)=+massFilteringScales(1)                                       &
            &                 *self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime)
       massFilteringScales(3)=+massFiltering
    end if
    return
  end subroutine gnedin2000ConditionsInitialODEs

  function gnedin2000ODEs(cosmologyParameters_,cosmologyFunctions_,linearGrowth_,time,massParticleMean,temperature,massFilteringODEs) result (massFilteringODEsRateOfChange)
    !!{RST
    Compute the rates of change of the filtering mass ODE system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigayear          , megaparsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    implicit none
    double precision                                         , dimension(3) :: massFilteringODEsRateOfChange
    class           (cosmologyParametersClass), intent(inout)               :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(inout)               :: cosmologyFunctions_
    class           (linearGrowthClass       ), intent(inout)               :: linearGrowth_
    double precision                          , intent(in   )               :: time                           , massParticleMean   , &
         &                                                                     temperature
    double precision                          , intent(in   ), dimension(3) :: massFilteringODEs
    double precision                                                        :: darkMatterFraction             , wavenumberFiltering, &
         &                                                                     wavenumberFilteringRateOfChange

    ! Compute dark matter mass fraction.
    darkMatterFraction              =1.0d0-cosmologyParameters_%OmegaBaryon()/cosmologyParameters_%OmegaMatter()
    ! Evaluate filtering mass composite terms.    
    !! These represent the 2nd order ODE given in equation 11 of Naoz & Barkana (2007, MNRAS, 377, 667;
    !! http://adsabs.harvard.edu/abs/2007MNRAS.377..667N), plus an ODE describing the evolution of the filtering mass in terms of
    !! the evolution of the growth factor and the filtering wavenumber.
    massFilteringODEsRateOfChange(1)=massFilteringODEs(2)
    if (massParticleMean > 0.0d0) then
       massFilteringODEsRateOfChange(2)=-2.0d0                                                                            &
            &                           *cosmologyFunctions_%expansionRate(cosmologyFunctions_%expansionFactor (time))    &
            &                           *massFilteringODEs(2)                                                             &
            &                           +darkMatterFraction                                                               &
            &                           /                                  cosmologyFunctions_%expansionFactor (time)**2  &
            &                           *boltzmannsConstant                                                               &
            &                           *temperature                                                                      &
            &                           /massParticleMean                                                                 &
            &                           *linearGrowth_%value                                                   (time)     &
            &                           *(                                                                                &
            &                             +1.0d0                                                                          &
            &                             +gnedin2000rLSS(                                                                &
            &                                                              cosmologyParameters_%OmegaMatter    (    )   , &
            &                                                              cosmologyFunctions_ %expansionFactor(time)     &
            &                                            )                                                                &
            &                            )                                                                                &
            &                           *gigayear  **2                                                                    &
            &                           /megaparsec**2
    else
       massFilteringODEsRateOfChange(2)=0.0d0
    end if
    if (massFilteringODEs(1) > 0.0d0) then
       wavenumberFiltering             =+sqrt(                           &
            &                                 +linearGrowth_%value(time) &
            &                                 /massFilteringODEs  (1   ) &
            &                                )
       wavenumberFilteringRateOfChange =+0.5d0                                                                                                 &
            &                           /sqrt(                                                                                                 &
            &                                 +linearGrowth_%value                                                                     (time)  &
            &                                 /massFilteringODEs                                                                       (1   )  &
            &                                )                                                                                                 &
            &                           *(                                                                                                     &
            &                             +linearGrowth_      %logarithmicDerivativeExpansionFactor                                    (time)  &
            &                             *cosmologyFunctions_%expansionRate                       (cosmologyFunctions_%expansionFactor(time)) &
            &                             *linearGrowth_      %value                                                                   (time)  &
            &                             *massFilteringODEs                                                                           (1   )  &
            &                             -linearGrowth_      %value                                                                   (time)  &
            &                             *massFilteringODEs                                                                           (2   )  &
            &                            )                                                                                                     &
            &                           /massFilteringODEs(1)**2
       massFilteringODEsRateOfChange(3)=-4.0d0                                     &
            &                           *Pi                                    **4 &
            &                           *cosmologyParameters_%densityCritical()    &
            &                           *cosmologyParameters_%OmegaMatter    ()    &
            &                           /wavenumberFiltering                   **4 &
            &                           *wavenumberFilteringRateOfChange
    else
       massFilteringODEsRateOfChange(3)=+0.0d0
    end if
    return
  end function gnedin2000ODEs

  function gnedin2000CoefficientsEarlyEpoch(self,time) result (coefficients)
    !!{RST
    Return the coefficients of the fitting function for the filtering mass at early epochs from :cite:t:`naoz_formation_2007`. Checks for valid range of redshift and cosmology for the fit to be valid.
    !!}
    use :: Error, only : Warn
    implicit none
    class           (intergalacticMediumFilteringMassGnedin2000), intent(inout) :: self
    double precision                                            , dimension(4)  :: coefficients
    double precision                                            , intent(in   ) :: time
    logical                                                     , save          :: warnedOmega =.false., warnedRedshift =.false.
    double precision                                                            :: omegaMatter         , expansionFactor        , &
         &                                                                         redshift

    ! Extract matter density and redshift.
    omegaMatter    =self%cosmologyParameters_%OmegaMatter                (               )
    expansionFactor=self%cosmologyFunctions_ %expansionFactor            (time           )
    redshift       =self%cosmologyFunctions_ %redshiftFromExpansionFactor(expansionFactor)
    ! Validate input.
    if (omegaMatter < 0.25d0 .or. omegaMatter >   0.40d0) then
       if (.not.warnedOmega) then
          !$omp critical (intergalacticMediumFilteringMassGnedin2000WarnedOmega)
          if (.not.warnedOmega) then
             call Warn('gnedin2000CoefficientsEarlyEpoch: matter density outside validated range of fitting function; 0.25 ≤ Ωₘ ≤ 0.40')
             warnedOmega=.true.
          end if
          !$omp end critical (intergalacticMediumFilteringMassGnedin2000WarnedOmega)
       end if
    end if
    if (redshift    < 7.00d0 .or. redshift    > 150.00d0) then
       if (.not.warnedRedshift) then
          !$omp critical (intergalacticMediumFilteringMassGnedin2000WarnedRedshift)
          if (.not.warnedRedshift) then
             call Warn('gnedin2000CoefficientsEarlyEpoch: redshift outside validated range of fitting function; 7 ≤ z ≤ 150'           )
             warnedRedshift=.true.
          end if
          !$omp end critical (intergalacticMediumFilteringMassGnedin2000WarnedRedshift)
       end if
    end if
    ! Evaluate fitting function.
    coefficients(1)=-0.38d0*(omegaMatter**2)+ 0.41d0*omegaMatter- 0.16d0
    coefficients(2)=+3.30d0*(omegaMatter**2)- 3.38d0*omegaMatter+ 1.15d0
    coefficients(3)=-9.64d0*(omegaMatter**2)+ 9.75d0*omegaMatter- 2.37d0
    coefficients(4)=+9.80d0*(omegaMatter**2)-10.68d0*omegaMatter+11.60d0
    return
  end function gnedin2000CoefficientsEarlyEpoch

  double precision function gnedin2000rLSS(omegaMatter,expansionFactor)
    !!{RST
    Evaluate the :math:`r_\mathrm{LSS}` parameter of :cite:t:`naoz_formation_2007` using their fitting formula.
    !!}
    implicit none
    double precision, intent(in   ) :: omegaMatter     , expansionFactor
    double precision                :: rLSSCoefficient1, rLSSCoefficient2, rLSSCoefficient3

    rLSSCoefficient1=1.0d-4*(-1.99d0*(omegaMatter**2)+2.41d0*omegaMatter+0.21d0)
    rLSSCoefficient2=1.0d-3*(+6.37d0*(omegaMatter**2)-6.99d0*omegaMatter-1.76d0)
    rLSSCoefficient3=1.0d-2*(-1.83d0*(omegaMatter**2)+2.40d0*omegaMatter-0.54d0)
    ! Note that the coefficients for the different exponents of expansion factor are reversed from that given in Naoz &
    ! Barkana. Without this change the fit for rLSS does not work.
    gnedin2000rLSS  =+rLSSCoefficient1/expansionFactor**1.5d0 &
         &           +rLSSCoefficient2/expansionFactor        &
         &           +rLSSCoefficient3
    return
  end function gnedin2000rLSS
