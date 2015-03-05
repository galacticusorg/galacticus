!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of dark matter halo mass accretion histories using the \cite{correa_accretion_2015} algorithm.

  !# <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryCorrea2015">
  !#  <description>Dark matter halo mass accretion histories using the \cite{correa_accretion_2015} algorithm.</description>
  !# </darkMatterHaloMassAccretionHistory>

  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryCorrea2015
     !% A dark matter halo mass accretion historiy class using the \cite{correa_accretion_2015} algorithm.
     private
   contains
     procedure :: time => correa2015Time
  end type darkMatterHaloMassAccretionHistoryCorrea2015
  
  interface darkMatterHaloMassAccretionHistoryCorrea2015
     !% Constructors for the {\normalfont \ttfamily correa2015} dark matter halo mass accretion history class.
     module procedure correa2015Constructor
  end interface darkMatterHaloMassAccretionHistoryCorrea2015
  
contains

  function correa2015Constructor()
    !% Generic constructor for the {\normalfont \ttfamily correa2015} dark matter halo mass accretion history class.
    implicit none
    type(darkMatterHaloMassAccretionHistoryCorrea2015), target :: correa2015Constructor
    
    return
  end function correa2015Constructor
  
  double precision function correa2015Time(self,node,mass)
    !% Compute the time corresponding to {\normalfont \ttfamily mass} in the mass accretion history of {\normalfont \ttfamily
    !% thisNode} using the algorithm of \cite{correa_accretion_2015}.
    use Root_Finder
    use Cosmology_Functions
    use Dark_Matter_Halos_Correa2015
    implicit none
    class           (darkMatterHaloMassAccretionHistoryCorrea2015), intent(inout)          :: self
    type            (treeNode                                    ), intent(inout), pointer :: node
    double precision                                              , intent(in   )          :: mass
    class           (nodeComponentBasic                          )               , pointer :: baseBasicComponent
    class           (cosmologyFunctionsClass                     )               , pointer :: cosmologyFunctions_
    double precision                                              , parameter              :: toleranceRelative  =1.0d-6
    double precision                                              , parameter              :: toleranceAbsolute  =0.0d0
    type            (rootFinder                                  ), save                   :: finder
    !$omp threadprivate(finder)
    double precision                                                                       :: baseRedshift               , baseTime, &
         &                                                                                    baseExpansionFactor        , baseMass, &
         &                                                                                    redshift                   , aTilde  , &
         &                                                                                    bTilde

    ! Get properties of the base node.
    baseBasicComponent => node%basic()
    baseMass=baseBasicComponent%mass()
    baseTime=baseBasicComponent%time()
    ! Get required objects.
    cosmologyFunctions_ => cosmologyFunctions()
    ! Determine the base redshift. 
    baseExpansionFactor=cosmologyFunctions_%expansionFactor            (baseTime           )
    baseRedshift       =cosmologyFunctions_%redshiftFromExpansionFactor(baseExpansionFactor)
    ! Find the a~ and b~ parameters.
    call Dark_Matter_Halo_Correa2015_Fit_Parameters(baseMass,baseExpansionFactor,aTilde,bTilde)
    ! Solve for the redshift corresponding to this mass.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(redshiftMassSolver                  )
       call finder%tolerance   (toleranceAbsolute ,toleranceRelative)
       call finder%rangeExpand (                                                           &
            &                   rangeExpandUpward          =2.0d0                        , &
            &                   rangeExpandType            =rangeExpandMultiplicative    , &
            &                   rangeExpandUpwardSignExpect=rangeExpandSignExpectNegative  &
            &)
    end if
    redshift=finder%find(rootRange=[baseRedshift,2.0d0/baseExpansionFactor-1.0d0])
    ! Convert redshift to time.
    correa2015Time=cosmologyFunctions_ %cosmicTime                 (           &
         &          cosmologyFunctions_%expansionFactorFromRedshift (          &
         &                                                            redshift &
         &                                                          )          &
         &                                                         )
    return

  contains

    double precision function redshiftMassSolver(redshift)
      !% Root solver function used in finding the redshift corresponding to a given mass in the \cite{correa_accretion_2015} mass
      !% accretion history algorithm.
      implicit none
      double precision, intent(in   ) :: redshift
      
      redshiftMassSolver=+baseMass                                          &
           &             *           (1.0d0+redshift-baseRedshift)**aTilde  &
           &             *exp(bTilde*(     +redshift-baseRedshift)        ) &
           &             -mass
      return
    end function redshiftMassSolver
    
  end function correa2015Time
