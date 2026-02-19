!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An implementation of dark matter halo mass accretion histories computed using the fitting function of \cite{diemer_splashback_2020}.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryDiemer2020">
   <description>Dark matter halo mass accretion histories computed using the fitting function of \cite{diemer_splashback_2020}.</description>
  </darkMatterHaloMassAccretionHistory>
  !!]
  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryDiemer2020
     !!{
     A dark matter halo mass accretion history class computed using the fitting function of \cite{diemer_splashback_2020}.
     !!}
     private
     class(cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
   contains
     final     ::                      diemer2020Destructor
     procedure :: time              => diemer2020Time
     procedure :: massAccretionRate => diemer2020MassAccretionRate
  end type darkMatterHaloMassAccretionHistoryDiemer2020

  interface darkMatterHaloMassAccretionHistoryDiemer2020
     !!{
     Constructors for the \refClass{darkMatterHaloMassAccretionHistoryDiemer2020} dark matter halo mass accretion history class.
     !!}
     module procedure diemer2020ConstructorParameters
     module procedure diemer2020ConstructorInternal
  end interface darkMatterHaloMassAccretionHistoryDiemer2020

contains

  function diemer2020ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloMassAccretionHistoryDiemer2020} dark matter halo mass accretion history class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterHaloMassAccretionHistoryDiemer2020)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                     ), pointer       :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass               ), pointer       :: cosmologicalMassVariance_
    class(criticalOverdensityClass                    ), pointer       :: criticalOverdensity_

    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=darkMatterHaloMassAccretionHistoryDiemer2020(cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function diemer2020ConstructorParameters

  function diemer2020ConstructorInternal(cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloMassAccretionHistoryDiemer2020} dark matter halo mass accretion history class.
    !!}
    implicit none
    type (darkMatterHaloMassAccretionHistoryDiemer2020)                        :: self
    class(cosmologyFunctionsClass                     ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass               ), intent(in   ), target :: cosmologicalMassVariance_
    class(criticalOverdensityClass                    ), intent(in   ), target :: criticalOverdensity_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function diemer2020ConstructorInternal

  subroutine diemer2020Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloMassAccretionHistoryDiemer2020} dark matter halo mass accretion history class.
    !!}
    implicit none
    type(darkMatterHaloMassAccretionHistoryDiemer2020), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%criticalOverdensity_"     />
    !!]
    return
  end subroutine diemer2020Destructor

  double precision function diemer2020Time(self,node,mass)
    !!{
    Compute the time corresponding to {\normalfont \ttfamily mass} in the mass accretion history.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (darkMatterHaloMassAccretionHistoryDiemer2020), intent(inout), target :: self
    type            (treeNode                                    ), intent(inout), target :: node
    double precision                                              , intent(in   )         :: mass
    !$GLC attributes unused :: self, node, mass

    diemer2020Time=0.0d0
    call Error_Report('"time" method is not supported'//{introspection:location})
    return
  end function diemer2020Time

  double precision function diemer2020Mass(self,node,time)
    !!{
    Compute the mass corresponding to {\normalfont \ttfamily time} in the mass accretion history.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (darkMatterHaloMassAccretionHistoryDiemer2020), intent(inout), target :: self
    type            (treeNode                                    ), intent(inout), target :: node
    double precision                                              , intent(in   )         :: time
    !$GLC attributes unused :: self, node, time

    diemer2020Mass=0.0d0
    call Error_Report('"mass" method is not supported'//{introspection:location})
    return
  end function diemer2020Mass

  double precision function diemer2020MassAccretionRate(self,node,time)
    !!{
    Compute the mass accretion rate at the given {\normalfont \ttfamily time} in the mass accretion history of
    {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloMassAccretionHistoryDiemer2020), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(in   ) :: time
    class           (nodeComponentBasic                          ), pointer       :: basic
    double precision                                                              :: peakHeight, expansionFactor      , &
         &                                                                           redshift  , massAccretionExponent, &
         &                                                                           A         , B

    basic => node%basic()
    ! Compute the peak height for this halo.
    peakHeight=+self%criticalOverdensity_     %value       (time=basic%time(),mass=basic%mass(),node=node) &
         &     /self%cosmologicalMassVariance_%rootVariance(time=basic%time(),mass=basic%mass()          )
    ! Get the expansion factor and redshift.
    expansionFactor=self%cosmologyFunctions_%expansionFactor            (basic%time           ())
    redshift       =self%cosmologyFunctions_%redshiftFromExpansionFactor(      expansionFactor  )
    ! Computing fitting function coefficients (equations 8 of Diemer 2020).
    A=+1.1721d0+0.3255d0*redshift
    B=-0.2565d0+0.0932d0*redshift-0.0571d0*redshift**2+0.0042d0*redshift**3
    ! Compute the mass accretion exponent (equation 7 of Diemer 2020).
    massAccretionExponent=+A*peakHeight        &
         &                +B*peakHeight**1.5d0
    ! Compute the mass accretion rate.
    diemer2020MassAccretionRate=+                          massAccretionExponent                  &
         &                      *basic                    %mass                 (               ) &
         &                      *self %cosmologyFunctions_%expansionRate        (expansionFactor)
    return
  end function diemer2020MassAccretionRate

