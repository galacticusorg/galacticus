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
  An implementation of dark matter halo mass accretion histories using the \cite{correa_accretion_2015} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Linear_Growth             , only : linearGrowthClass
  use :: Root_Finder               , only : rootFinder

  !![
  <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryCorrea2015">
   <description>Dark matter halo mass accretion histories using the \cite{correa_accretion_2015} algorithm.</description>
  </darkMatterHaloMassAccretionHistory>
  !!]
  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryCorrea2015
     !!{
     A dark matter halo mass accretion history class using the \cite{correa_accretion_2015} algorithm.
     !!}
     private
     class(cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class(linearGrowthClass            ), pointer :: linearGrowth_             => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     type (rootFinder                   )          :: finder
   contains
     final     ::                      correa2015Destructor
     procedure :: time              => correa2015Time
     procedure :: massAccretionRate => correa2015MassAccretionRate
  end type darkMatterHaloMassAccretionHistoryCorrea2015

  interface darkMatterHaloMassAccretionHistoryCorrea2015
     !!{
     Constructors for the \refClass{darkMatterHaloMassAccretionHistoryCorrea2015} dark matter halo mass accretion history class.
     !!}
     module procedure correa2015ConstructorParameters
     module procedure correa2015ConstructorInternal
  end interface darkMatterHaloMassAccretionHistoryCorrea2015

contains

  function correa2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloMassAccretionHistoryCorrea2015} dark matter halo mass accretion history class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterHaloMassAccretionHistoryCorrea2015)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                     ), pointer       :: cosmologyFunctions_
    class(linearGrowthClass                           ), pointer       :: linearGrowth_
    class(cosmologicalMassVarianceClass               ), pointer       :: cosmologicalMassVariance_

    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=darkMatterHaloMassAccretionHistoryCorrea2015(cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function correa2015ConstructorParameters

  function correa2015ConstructorInternal(cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloMassAccretionHistoryCorrea2015} dark matter halo mass accretion history class.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative
    implicit none
    type            (darkMatterHaloMassAccretionHistoryCorrea2015)                        :: self
    class           (cosmologyFunctionsClass                     ), intent(in   ), target :: cosmologyFunctions_
    class           (linearGrowthClass                           ), intent(in   ), target :: linearGrowth_
    class           (cosmologicalMassVarianceClass               ), intent(in   ), target :: cosmologicalMassVariance_
    double precision                                              , parameter             :: toleranceRelative        =1.0d-6,toleranceAbsolute=0.0d+0
    !![
    <constructorAssign variables="*cosmologyFunctions_, *linearGrowth_, *cosmologicalMassVariance_"/>
    !!]
    
    self%finder=rootFinder(                                                           &
         &                 toleranceAbsolute          =toleranceAbsolute            , &
         &                 toleranceRelative          =toleranceRelative            , &
         &                 rangeExpandUpward          =2.0d0                        , &
         &                 rangeExpandType            =rangeExpandMultiplicative    , &
         &                 rangeExpandUpwardSignExpect=rangeExpandSignExpectNegative  &
         &                )
    return
  end function correa2015ConstructorInternal

  subroutine correa2015Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloMassAccretionHistoryCorrea2015} dark matter halo mass accretion history class.
    !!}
    implicit none
    type(darkMatterHaloMassAccretionHistoryCorrea2015), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine correa2015Destructor

  double precision function correa2015Time(self,node,mass)
    !!{
    Compute the time corresponding to {\normalfont \ttfamily mass} in the mass accretion history of {\normalfont \ttfamily
    node} using the algorithm of \cite{correa_accretion_2015}.
    !!}
    use :: Dark_Matter_Halos_Correa2015, only : Dark_Matter_Halo_Correa2015_Fit_Parameters
    use :: Galacticus_Nodes            , only : nodeComponentBasic                        , treeNode
    implicit none
    class           (darkMatterHaloMassAccretionHistoryCorrea2015), intent(inout), target :: self
    type            (treeNode                                    ), intent(inout), target :: node
    double precision                                              , intent(in   )         :: mass
    class           (nodeComponentBasic                          ), pointer               :: basicBase
    double precision                                                                      :: baseRedshift               , baseTime, &
         &                                                                                   baseExpansionFactor        , baseMass, &
         &                                                                                   redshift                   , aTilde  , &
         &                                                                                   bTilde
    !$GLC attributes unused :: self

    ! Get properties of the base node.
    basicBase => node     %basic()
    baseMass  =  basicBase%mass ()
    baseTime  =  basicBase%time ()
    ! Determine the base redshift.
    baseExpansionFactor=self%cosmologyFunctions_%expansionFactor            (baseTime           )
    baseRedshift       =self%cosmologyFunctions_%redshiftFromExpansionFactor(baseExpansionFactor)
    ! Find the a~ and b~ parameters.
    call Dark_Matter_Halo_Correa2015_Fit_Parameters(baseMass,baseExpansionFactor,self%cosmologyFunctions_,self%linearGrowth_,self%cosmologicalMassVariance_,aTilde,bTilde)
    ! Solve for the redshift corresponding to this mass.
    call self%finder%rootFunction(redshiftMassSolver)
    redshift=self%finder%find(rootRange=[baseRedshift,2.0d0/baseExpansionFactor-1.0d0])
    ! Convert redshift to time.
    correa2015Time=self%cosmologyFunctions_%cosmicTime                 (          &
         &         self%cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                               redshift &
         &                                                              )         &
         &                                                             )
    return

  contains

    double precision function redshiftMassSolver(redshift)
      !!{
      Root solver function used in finding the redshift corresponding to a given mass in the \cite{correa_accretion_2015} mass
      accretion history algorithm.
      !!}
      implicit none
      double precision, intent(in   ) :: redshift

      redshiftMassSolver=+baseMass                                          &
           &             *           (1.0d0+redshift-baseRedshift)**aTilde  &
           &             *exp(bTilde*(     +redshift-baseRedshift)        ) &
           &             -mass
      return
    end function redshiftMassSolver

  end function correa2015Time

  double precision function correa2015MassAccretionRate(self,node,time)
    !!{
    Compute the mass accretion rate at the given time {\normalfont \ttfamily mass} in the mass accretion history of
    {\normalfont \ttfamily node} using the algorithm of \cite{correa_accretion_2015}.
    !!}
    use :: Dark_Matter_Halos_Correa2015, only : Dark_Matter_Halo_Correa2015_Fit_Parameters
    use :: Galacticus_Nodes            , only : nodeComponentBasic                        , treeNode
    implicit none
    class           (darkMatterHaloMassAccretionHistoryCorrea2015), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(in   ) :: time
    class           (nodeComponentBasic                          ), pointer       :: basicBase
    double precision                                                              :: baseRedshift       , baseTime                 , &
         &                                                                           baseExpansionFactor, baseMass                 , &
         &                                                                           redshift           , aTilde                   , &
         &                                                                           bTilde             , expansionFactor          , &
         &                                                                           redshift           , massAccretionRateRedshift
    !$GLC attributes unused :: self

    ! Get properties of the base node.
    basicBase => node     %basic()
    baseMass  =  basicBase%mass ()
    baseTime  =  basicBase%time ()
    ! Determine the base redshift.
    baseExpansionFactor=self%cosmologyFunctions_%expansionFactor            (baseTime           )
    baseRedshift       =self%cosmologyFunctions_%redshiftFromExpansionFactor(baseExpansionFactor)
    ! Find the a~ and b~ parameters.
    call Dark_Matter_Halo_Correa2015_Fit_Parameters(baseMass,baseExpansionFactor,self%cosmologyFunctions_,self%linearGrowth_,self%cosmologicalMassVariance_,aTilde,bTilde)
    ! Compute the mass accretion rate per unit redshift.
    expansionFactor          =+self%cosmologyFunctions_%expansionFactor            (time           )
    redshift                 =+self%cosmologyFunctions_%redshiftFromExpansionFactor(expansionFactor)
    massAccretionRateRedshift=+baseMass                                           &
         &                    *            (1.0d0+redshift-baseRedshift) **aTilde &
         &                    *exp(+bTilde*(     +redshift-baseRedshift))         &
         &                    *   (                                               &
         &                         +aTilde/(1.0d0+redshift-baseRedshift)          &
         &                         +bTilde                                        &
         &                        )
    ! Convert to mass accretion rate per unit time.
    correa2015MassAccretionRate=-massAccretionRateRedshift                               &
         &                      *self%cosmologyFunctions_%expansionRate(expansionFactor) &
         &                      /                                       expansionFactor
    return
  end function correa2015MassAccretionRate
  
