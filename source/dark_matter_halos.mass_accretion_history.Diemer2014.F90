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

  !% An implementation of dark matter halo mass accretion histories computed using the fitting function of \cite{diemer_dependence_2014}.

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !# <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryDiemer2014">
  !#  <description>Dark matter halo mass accretion histories computed using the fitting function of \cite{diemer_dependence_2014}.</description>
  !# </darkMatterHaloMassAccretionHistory>
  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryDiemer2014
     !% A dark matter halo mass accretion history class computed using the fitting function of \cite{diemer_dependence_2014}.
     private
     class(cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
   contains
     final     ::                      diemer2014Destructor
     procedure :: time              => diemer2014Time
     procedure :: massAccretionRate => diemer2014MassAccretionRate
  end type darkMatterHaloMassAccretionHistoryDiemer2014

  interface darkMatterHaloMassAccretionHistoryDiemer2014
     !% Constructors for the {\normalfont \ttfamily diemer2014} dark matter halo mass accretion history class.
     module procedure diemer2014ConstructorParameters
     module procedure diemer2014ConstructorInternal
  end interface darkMatterHaloMassAccretionHistoryDiemer2014

contains

  function diemer2014ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily diemer2014} dark matter halo mass accretion history class which takes a
    !% parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterHaloMassAccretionHistoryDiemer2014)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                     ), pointer       :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass               ), pointer       :: cosmologicalMassVariance_
    class(criticalOverdensityClass                    ), pointer       :: criticalOverdensity_

    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    self=darkMatterHaloMassAccretionHistoryDiemer2014(cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"      />
    !# <objectDestructor name="criticalOverdensity_"     />
    !# <objectDestructor name="cosmologicalMassVariance_"/>
    return
  end function diemer2014ConstructorParameters

  function diemer2014ConstructorInternal(cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !% Generic constructor for the {\normalfont \ttfamily diemer2014} dark matter halo mass accretion history class.
    implicit none
    type (darkMatterHaloMassAccretionHistoryDiemer2014)                        :: self
    class(cosmologyFunctionsClass                     ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass               ), intent(in   ), target :: cosmologicalMassVariance_
    class(criticalOverdensityClass                    ), intent(in   ), target :: criticalOverdensity_
    !# <constructorAssign variables="*cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_"/>

    return
  end function diemer2014ConstructorInternal

  subroutine diemer2014Destructor(self)
    !% Destructor for the {\normalfont \ttfamily diemer2014} dark matter halo mass accretion history class.
    implicit none
    type(darkMatterHaloMassAccretionHistoryDiemer2014), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"      />
    !# <objectDestructor name="self%cosmologicalMassVariance_"/>
    !# <objectDestructor name="self%criticalOverdensity_"     />
    return
  end subroutine diemer2014Destructor

  double precision function diemer2014Time(self,node,mass)
    !% Compute the time corresponding to {\normalfont \ttfamily mass} in the mass accretion history.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (darkMatterHaloMassAccretionHistoryDiemer2014), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(in   ) :: mass
    !$GLC attributes unused :: self, node, mass

    diemer2014Time=0.0d0
    call Galacticus_Error_Report('"time" method is not supported'//{introspection:location})
    return
  end function diemer2014Time

  double precision function diemer2014MassAccretionRate(self,node,time)
    !% Compute the mass accretion rate at the given {\normalfont \ttfamily time} in the mass accretion history of
    !% {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloMassAccretionHistoryDiemer2014), intent(inout) :: self
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
    ! Computing fitting function coefficients (equations 8 of Diemer & Kravtsov 2014).
    A=+1.1721d0+0.3255*redshift
    B=-0.2565d0+0.0932*redshift-0.0571*redshift**2+0.0042*redshift**3
    ! Compute the mass accretion exponent (equation 7 of Diemer & Kravtsov 2014).
    massAccretionExponent=+A*peakHeight        &
         &                +B*peakHeight**1.5d0
    ! Compute the mass accretion rate.
    diemer2014MassAccretionRate=+                          massAccretionExponent                  &
         &                      *basic                    %mass                 (               ) &
         &                      *self %cosmologyFunctions_%expansionRate        (expansionFactor)
    return
  end function diemer2014MassAccretionRate

