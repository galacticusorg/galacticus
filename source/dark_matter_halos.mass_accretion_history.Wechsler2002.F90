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
  An implementation of dark matter halo mass accretion histories using the \cite{wechsler_concentrations_2002} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Kind_Numbers              , only : kind_int8

  !![
  <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryWechsler2002">
   <description>
    A dark matter halo mass accretion history class in which the mass accretion history is given by
    \citep{wechsler_concentrations_2002}:
    \begin{equation}
    M(t) = M(t_0) \exp \left( - 2 a_\mathrm{c} \left[ {a(t_0)\over a(t)}-1 \right] \right),
    \end{equation}
    where $t_0$ is some reference time and $a_\mathrm{c}$ is a characteristic expansion factor defined by
    \cite{wechsler_concentrations_2002} to correspond to the formation time of the halo (using the formation time definition of
    \citealt{bullock_profiles_2001}).
   </description>
  </darkMatterHaloMassAccretionHistory>
  !!]
  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryWechsler2002
     !!{
     A dark matter halo mass accretion history class using the \cite{wechsler_concentrations_2002} algorithm.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     logical                                                  :: formationRedshiftCompute
     double precision                                         :: formationRedshift
     integer         (kind=kind_int8               )          :: lastUniqueID              =  -huge(0_kind_int8)
     double precision                                         :: timeFormationPrevious     =  -huge(0.0d0      ), massPrevious=-huge(0.0d0)
   contains
     !![
     <methods>
       <method description="Reset memoized calculations."            method="calculationReset"          />
       <method description="Compute the formation expansion factor." method="expansionFactorAtFormation"/>
     </methods>
     !!]
     final     ::                               wechsler2002Destructor
     procedure :: autoHook                   => wechsler2002AutoHook
     procedure :: calculationReset           => wechsler2002CalculationReset
     procedure :: time                       => wechsler2002Time
     procedure :: massAccretionRate          => wechsler2002MassAccretionRate
     procedure :: expansionFactorAtFormation => wechsler2002ExpansionFactorAtFormation
  end type darkMatterHaloMassAccretionHistoryWechsler2002

  interface darkMatterHaloMassAccretionHistoryWechsler2002
     !!{
     Constructors for the \refClass{darkMatterHaloMassAccretionHistoryWechsler2002} dark matter halo mass accretion history class.
     !!}
     module procedure wechsler2002ConstructorParameters
     module procedure wechsler2002ConstructorInternal
  end interface darkMatterHaloMassAccretionHistoryWechsler2002

contains

  function wechsler2002ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily wechsler2002} dark matter halo mass accretion history class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterHaloMassAccretionHistoryWechsler2002)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass                      ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                 ), pointer       :: cosmologicalMassVariance_
    logical                                                                         :: formationRedshiftCompute
    double precision                                                                :: formationRedshift

    !![
    <inputParameter>
      <name>formationRedshiftCompute</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, compute formation redshift automatically for \cite{wechsler_concentrations_2002} halo mass accretion histories.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (.not.formationRedshiftCompute) then
       ! In this case, read the formation redshift.
       !![
       <inputParameter>
         <name>formationRedshift</name>
         <description>The formation redshift to use in \cite{wechsler_concentrations_2002} halo mass accretion histories.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <conditionalCall>
     <call>self=darkMatterHaloMassAccretionHistoryWechsler2002(cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,formationRedshiftCompute{conditions})</call>
     <argument name="formationRedshift" value="formationRedshift" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function wechsler2002ConstructorParameters

  function wechsler2002ConstructorInternal(cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,formationRedshiftCompute,formationRedshift) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloMassAccretionHistoryWechsler2002} dark matter halo mass accretion history class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (darkMatterHaloMassAccretionHistoryWechsler2002)                          :: self
    class           (cosmologyFunctionsClass                       ), intent(in   ), target   :: cosmologyFunctions_
    class           (criticalOverdensityClass                      ), intent(in   ), target   :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                 ), intent(in   ), target   :: cosmologicalMassVariance_
    logical                                                         , intent(in   )           :: formationRedshiftCompute
    double precision                                                , intent(in   ), optional :: formationRedshift
    !![
    <constructorAssign variables="*cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, formationRedshiftCompute, formationRedshift"/>
    !!]

    if (.not.formationRedshiftCompute.and..not.present(formationRedshift)) call Error_Report('formation redshift must be provided'//{introspection:location})
    return
  end function wechsler2002ConstructorInternal

  subroutine wechsler2002AutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterHaloMassAccretionHistoryWechsler2002), intent(inout) :: self

    call calculationResetEvent%attach(self,wechsler2002CalculationReset,openMPThreadBindingAllLevels,label='darkMatterHaloMassAccretionHistoryWechsler2002')
    return
  end subroutine wechsler2002AutoHook
  
  subroutine wechsler2002Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloMassAccretionHistoryWechsler2002} dark matter halo mass accretion history class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterHaloMassAccretionHistoryWechsler2002), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,wechsler2002CalculationReset)) call calculationResetEvent%detach(self,wechsler2002CalculationReset)
    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine wechsler2002Destructor

  subroutine wechsler2002CalculationReset(self,node,uniqueID)
    !!{
    Reset the cooling radius calculation.
    !!}
    implicit none
    class  (darkMatterHaloMassAccretionHistoryWechsler2002), intent(inout) :: self
    type   (treeNode                                      ), intent(inout) :: node
    integer(kind_int8                                     ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%timeFormationPrevious=-huge(0.0d0)
    self%massPrevious         =-huge(0.0d0)
    self%lastUniqueID         =uniqueID
    return
  end subroutine wechsler2002CalculationReset

  double precision function wechsler2002Time(self,node,mass)
    !!{
    Compute the time corresponding to {\normalfont \ttfamily mass} in the mass accretion history of {\normalfont \ttfamily node} using the algorithm of
    \cite{wechsler_concentrations_2002}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterHaloMassAccretionHistoryWechsler2002), intent(inout), target :: self
    type            (treeNode                                      ), intent(inout), target :: node
    double precision                                                , intent(in   )         :: mass
    class           (nodeComponentBasic                            ), pointer               :: basicBase
    double precision                                                                        :: expansionFactor                   , expansionFactorBase, &
         &                                                                                     mergerTreeFormationExpansionFactor

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Compute formation time if necessary.
    if (self%massPrevious /= mass) then
       basicBase => node%basic()
       select case (self%formationRedshiftCompute)
       case (.true.)
          ! Compute the expansion factor at formation.
          mergerTreeFormationExpansionFactor=self%expansionFactorAtFormation(basicBase%mass(),node)
       case (.false.)
          ! Use the specified formation redshift.
          mergerTreeFormationExpansionFactor=self%cosmologyFunctions_%expansionFactorFromRedshift(self%formationRedshift)
       end select
       ! Get the expansion factor at the tree base.
       expansionFactorBase=self%cosmologyFunctions_%expansionFactor(basicBase%time())
       ! Compute the expansion factor for the current node.
       expansionFactor=expansionFactorBase/(1.0d0-0.5d0*log(mass/basicBase%mass())/mergerTreeFormationExpansionFactor)
       ! Find the time corresponding to this expansion factor.
       self%timeFormationPrevious=self%cosmologyFunctions_%cosmicTime(expansionFactor)
       self%massPrevious         =mass
    end if
    wechsler2002Time=self%timeFormationPrevious
    return
  end function wechsler2002Time

  double precision function wechsler2002MassAccretionRate(self,node,time)
    !!{
    Compute the mass accretion rate at the given {\normalfont \ttfamily time} in the mass accretion history of {\normalfont
    \ttfamily node} using the algorithm of \cite{wechsler_concentrations_2002}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterHaloMassAccretionHistoryWechsler2002), intent(inout) :: self
    type            (treeNode                                      ), intent(inout) :: node
    double precision                                                , intent(in   ) :: time
    class           (nodeComponentBasic                            ), pointer       :: basicBase
    double precision                                                                :: expansionFactor                   , expansionFactorBase, &
         &                                                                             mergerTreeFormationExpansionFactor

    basicBase => node%basic()
    select case (self%formationRedshiftCompute)
    case (.true.)
       ! Compute the expansion factor at formation.
       mergerTreeFormationExpansionFactor=self%expansionFactorAtFormation(basicBase%mass(),node)
       if (mergerTreeFormationExpansionFactor == huge(0.0d0)) then
          ! No collapse time was found - the accretion rate is indeterminate.
          wechsler2002MassAccretionRate=-huge(0.0d0)
          return
       end if
    case (.false.)
       ! Use the specified formation redshift.
       mergerTreeFormationExpansionFactor=self%cosmologyFunctions_%expansionFactorFromRedshift(self%formationRedshift)
    end select
    ! Get the expansion factor at the tree base.
    expansionFactorBase=self%cosmologyFunctions_%expansionFactor(basicBase%time())
    ! Get the expansion factor at the current time.
    expansionFactor    =self%cosmologyFunctions_%expansionFactor(          time  )
    ! Compute the mass accretion rate.
    wechsler2002MassAccretionRate=+basicBase%mass()                                        &
         &                        *2.0d0                                                   &
         &                        *mergerTreeFormationExpansionFactor                      &
         &                        *expansionFactorBase                                     &
         &                        /expansionFactor                                         &
         &                        *exp(                                                    &
         &                             -2.0d0                                              &
         &                             *mergerTreeFormationExpansionFactor                 &
         &                             *(                                                  &
         &                               +expansionFactorBase                              &
         &                               /expansionFactor                                  &
         &                               -1.0d0                                            &
         &                              )                                                  &
         &                            )                                                    &
         &                        *self%cosmologyFunctions_%expansionRate(expansionFactor)
    return
  end function wechsler2002MassAccretionRate

  double precision function wechsler2002ExpansionFactorAtFormation(self,haloMass,node)
    !!{
    Computes the expansion factor at formation using the simple model of \cite{bullock_profiles_2001}.
    !!}
    use :: Error, only : Error_Report, errorStatusSuccess, errorStatusOutOfRange
    implicit none
    class           (darkMatterHaloMassAccretionHistoryWechsler2002), intent(inout) :: self
    type            (treeNode                                      ), intent(inout) :: node
    double precision                                                , intent(in   ) :: haloMass
    double precision                                                , parameter     :: haloMassFraction   =0.015d0 ! Wechsler et al. (2002;  Astrophysical Journal, 568:52-70).
    double precision                                                                :: formationTime              , haloMassCharacteristic, &
         &                                                                             sigmaCharacteristic
    integer                                                                         :: status

    ! Compute the characteristic mass at formation time.
    haloMassCharacteristic=haloMassFraction*haloMass
    ! Compute the corresponding rms fluctuation in the density field (i.e. σ(M)).
    sigmaCharacteristic=self%cosmologicalMassVariance_%rootVariance(haloMassCharacteristic,self%cosmologyFunctions_%cosmicTime(1.0d0))
    ! Get the time at which this equals the critical overdensity for collapse.
    formationTime=self%criticalOverdensity_%timeOfCollapse(criticalOverdensity=sigmaCharacteristic,mass=haloMass,node=node,status=status)
    select case (status)
    case (errorStatusSuccess)
       ! Get the corresponding expansion factor.
       wechsler2002ExpansionFactorAtFormation=self%cosmologyFunctions_%expansionFactor(formationTime)
    case (errorStatusOutOfRange)
       ! Collapse never occurs for this halo mass.
       wechsler2002ExpansionFactorAtFormation=huge(0.0d0)
    case default
       ! Unknown error.
       wechsler2002ExpansionFactorAtFormation=huge(0.0d0)
       call Error_Report('failed to find expansion factor at formation'//{introspection:location})
    end select
    return
  end function wechsler2002ExpansionFactorAtFormation
