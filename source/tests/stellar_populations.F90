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
Contains a program to test stellar populations.
!!}

program Test_Stellar_Populations
  !!{RST
  Tests of stellar populations.
  !!}
  use :: Abundances_Structure                      , only : abundances
  use :: Display                                   , only : displayVerbositySet                    , verbosityLevelWorking
  use :: Events_Hooks                              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities                , only : Functions_Global_Set
  use :: Input_Paths                               , only : inputPath                              , pathTypeDataStatic
  use :: ISO_Varying_String                        , only : char                                   , varying_string         , assignment(=)       , operator(//)
  use :: Input_Parameters                          , only : inputParameters
  use :: Galacticus_Nodes                          , only : nodeClassHierarchyInitialize
  use :: Node_Components                           , only : Node_Components_Initialize
  use :: Numerical_Constants_Astronomical          , only : metallicitySolar
  use :: Numerical_Integration2                    , only : integratorCompositeGaussKronrod1D
  use :: Stellar_Astrophysics                      , only : stellarAstrophysics                    , stellarAstrophysicsFile
  use :: Stellar_Astrophysics_Tracks               , only : stellarTracksFile
  use :: Stellar_Astrophysics_Winds                , only : stellarWindsLeitherer1992
  use :: Stellar_Feedback                          , only : stellarFeedbackStandard
  use :: Stellar_Population_Spectra                , only : stellarPopulationSpectraFSPS
  use :: Stellar_Populations                       , only : stellarPopulationStandard
  use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionChabrier2001
  use :: Supernovae_Population_III                 , only : supernovaePopulationIIIHegerWoosley2002
  use :: Supernovae_Type_Ia                        , only : supernovaeTypeIaNagashima2005
  use :: Unit_Tests                                , only : Assert                                 , Unit_Tests_Begin_Group , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                                         , parameter :: ageMinimum               =0.0d0, ageMaximum=10.0d0
  type            (inputParameters                        ), target    :: parameters
  type            (abundances                             )            :: abundances_
  type            (initialMassFunctionChabrier2001        )            :: initialMassFunction_
  type            (stellarAstrophysicsFile                )            :: stellarAstrophysics_
  type            (stellarPopulationStandard              )            :: stellarPopulation_
  type            (stellarFeedbackStandard                )            :: stellarFeedback_
  type            (stellarTracksFile                      )            :: stellarTracks_
  type            (stellarWindsLeitherer1992              )            :: stellarWinds_
  type            (supernovaeTypeIaNagashima2005          )            :: supernovaeTypeIa_
  type            (supernovaePopulationIIIHegerWoosley2002)            :: supernovaePopulationIII_
  type            (stellarPopulationSpectraFSPS           )            :: stellarPopulationSpectra_
  double precision                                                     :: recycledMass                   , yieldMetals
  type            (varying_string                         )            :: fileNameCompilation            , fileNameTracks
  ! A point at which the tabulated recycled fraction and yield are evaluated without interpolation. Both coordinates are nodes
  ! of the age and metallicity grids on which the tables are built (50 points logarithmically spaced over 10⁻³ to 10² Gyr, and
  ! zero followed by 9 points logarithmically spaced over 10⁻⁴ to 6 × 10⁻²), so the linear interpolation in each axis returns
  ! the tabulated value exactly.
  double precision                                         , parameter :: ageNode                        =1.2067926406393298d+01
  double precision                                         , parameter :: metallicityNode                =2.6970086794142582d-02
  type            (abundances                             )            :: abundancesNode
  type            (integratorCompositeGaussKronrod1D       )            :: integrator_
  double precision                                                     :: recycledFractionTabulated      , recycledFractionIntegrated, &
       &                                                                  yieldTabulated                 , yieldIntegrated
  
  call displayVerbositySet(verbosityLevelWorking)
  parameters=inputParameters()
  call Functions_Global_Set        (          )
  call eventsHooksInitialize       (          )
  call nodeClassHierarchyInitialize(parameters)
  call Node_Components_Initialize  (parameters)
  call abundances_%metallicitySet  (metallicitySolar)
  fileNameCompilation      =inputPath(pathTypeDataStatic)//'stellarAstrophysics/stellarPropertiesCompilationStandard.xml'
  fileNameTracks           =inputPath(pathTypeDataStatic)//'stellarAstrophysics/Stellar_Tracks_Padova.hdf5'
  initialMassFunction_     =initialMassFunctionChabrier2001        (                                                                 &
       &                                                            massLower                            =+0.10d0                  , &
       &                                                            massTransition                       =+1.00d0                  , &
       &                                                            massUpper                            =+1.25d2                  , &
       &                                                            exponent                             =-2.30d0                  , &
       &                                                            massCharacteristic                   =+0.08d0                  , &
       &                                                            sigma                                =+0.69d0                    &
       &                                                           )
  stellarAstrophysics_     =stellarAstrophysicsFile                (                                                                 &
       &                                                            fileName                             =char(fileNameCompilation)  &
       &                                                           )
  stellarTracks_           =stellarTracksFile                      (                                                                 &
       &                                                            fileName                             =char(fileNameTracks     )  &
       &                                                           )
  stellarWinds_            =stellarWindsLeitherer1992              (                                                                 &
       &                                                            stellarTracks_                       =stellarTracks_             &
       &                                                           )
  supernovaeTypeIa_        =supernovaeTypeIaNagashima2005          (                                                                 &
       &                                                            stellarAstrophysics_                 =stellarAstrophysics_     , &
       &                                                            fileName                             =char(inputPath(pathTypeDataStatic))//'stellarAstrophysics/Supernovae_Type_Ia_Yields.xml'  &
       &                                                           )
  supernovaePopulationIII_ =supernovaePopulationIIIHegerWoosley2002(                                                                 &
       &                                                            stellarAstrophysics_                 =stellarAstrophysics_       &
       &                                                           )
  stellarFeedback_         =stellarFeedbackStandard                (                                                                 &
       &                                                            initialMassForSupernovaeTypeII       =8.0d00                   , &
       &                                                            supernovaEnergy                      =1.0d51                   , &
       &                                                            supernovaeTypeIa_                    =supernovaeTypeIa_        , &
       &                                                            supernovaePopulationIII_             =supernovaePopulationIII_ , &
       &                                                            stellarWinds_                        =stellarWinds_            , &
       &                                                            stellarAstrophysics_                 =stellarAstrophysics_       &
       &                                                           )
  stellarPopulationSpectra_=stellarPopulationSpectraFSPS           (                                                                 &
       &                                                            forceZeroMetallicity                 =.false.                  , &
       &                                                            initialMassFunction_                 =initialMassFunction_       &
       &                                                           )
  stellarPopulation_       =stellarPopulationStandard              (                                                                 &
       &                                                            instantaneousRecyclingApproximation  =.false.                  , &
       &                                                            instantaneousYieldApproximation      =.false.                  , &
       &                                                            instantaneousEnergyInputApproximation=.false.                  , &
       &                                                            massLongLived                        =1.0d0                    , &
       &                                                            ageEffective                         =1.0d1                    , &
       &                                                            initialMassFunction_                 =initialMassFunction_     , &
       &                                                            stellarAstrophysics_                 =stellarAstrophysics_     , &
       &                                                            stellarFeedback_                     =stellarFeedback_         , &
       &                                                            supernovaeTypeIa_                    =supernovaeTypeIa_        , &
       &                                                            stellarPopulationSpectra_            =stellarPopulationSpectra_  &
       &                                                           )
  call Unit_Tests_Begin_Group("Stellar population functions")
  ! Compute the recycled mass, yield etc. by getting the mean production rate in our time interval and multiplying by the size of
  ! the interval. This is compared to a previously computed value for the Chabrier (2001) IMF.
  recycledMass=stellarPopulation_%recycledFractionInstantaneous()
  call Assert('recycled fraction',recycledMass,4.6d-1,relTol=1.0d-2)
  yieldMetals =stellarPopulation_%yieldInstantaneous    ()
  call Assert('metal yield'      ,yieldMetals ,3.7d-2,relTol=1.0d-2)
  call Unit_Tests_End_Group  ()
  ! Recycled fraction and metal yield, integrated independently over the initial mass function and compared against the values
  ! the stellar population class returns.
  !
  ! The two assertions above compare against values recorded from a previous run of this code, so they test that the result has
  ! not changed, not that it is right. The assertions below instead assemble the same quantities from their definitions,
  !
  !   f_recycled = ∫ φ(M) M_ejected(M,Z) dM   over stars whose lifetime is less than the age,
  !   y_metals   = ∫ φ(M) M_yield  (M,Z) dM   over the same, plus the Type Ia contribution over *all* masses,
  !
  ! taking the ejected masses, yields and lifetimes from Galacticus itself. What is checked is therefore the assembly - the
  ! restriction to evolved stars, the weighting by the initial mass function, and the fact that the Type Ia term is deliberately
  ! *not* restricted to evolved stars - rather than the tabulated stellar data, which both sides share.
  !
  ! The comparison is made at a node of the age and metallicity grids, where the class performs no interpolation. This matters:
  ! the tables are coarse in metallicity, and at Solar metallicity - which falls 45% of the way through its interval - linear
  ! interpolation of the recycled fraction differs from a higher-order estimate by 2.2%, and of the yield by more. Comparing
  ! there would measure the metallicity grid rather than the physics. At a node the only difference expected is the quadrature
  ! itself, for which the tables are built to a relative tolerance of 10⁻⁴ (recycled fraction) and 10⁻⁵ (yields); the
  ! assertions below allow 10⁻³.
  call Unit_Tests_Begin_Group("Recycled fraction and yield, assembled independently")
  call abundancesNode%metallicitySet(metallicityNode)
  recycledFractionTabulated =stellarPopulation_%rateRecycling(abundancesNode,ageMinimum=0.0d0,ageMaximum=ageNode)*ageNode
  yieldTabulated            =stellarPopulation_%rateYield    (abundancesNode,ageMinimum=0.0d0,ageMaximum=ageNode)*ageNode
  call integrator_%initialize  (24    ,61    )
  call integrator_%toleranceSet(1.0d-3,1.0d-4)
  call integrator_%integrandSet(integrandRecycledFraction)
  recycledFractionIntegrated=integrator_%evaluate(initialMassFunction_%massMinimum(),initialMassFunction_%massMaximum())
  call integrator_%toleranceSet(1.0d-4,1.0d-5)
  call integrator_%integrandSet(integrandYield           )
  yieldIntegrated           =integrator_%evaluate(initialMassFunction_%massMinimum(),initialMassFunction_%massMaximum())
  call Assert('recycled fraction at a grid node',recycledFractionIntegrated,recycledFractionTabulated,relTol=1.0d-3)
  call Assert('metal yield at a grid node'      ,yieldIntegrated           ,yieldTabulated           ,relTol=1.0d-3)
  call Unit_Tests_End_Group  ()
  call Unit_Tests_Finish     ()

contains

  double precision function integrandRecycledFraction(massInitial)
    !!{RST
    Integrand giving the mass returned to the interstellar medium per unit mass formed, by stars which have evolved off of the
    main sequence by the age under consideration.
    !!}
    implicit none
    double precision, intent(in   ) :: massInitial

    if (stellarAstrophysics_%lifetime(massInitial,metallicityNode) < ageNode) then
       integrandRecycledFraction=+initialMassFunction_%phi        (massInitial                ) &
            &                    *stellarAstrophysics_%massEjected(massInitial,metallicityNode)
    else
       integrandRecycledFraction=+0.0d0
    end if
    return
  end function integrandRecycledFraction

  double precision function integrandYield(massInitial)
    !!{RST
    Integrand giving the mass of metals returned to the interstellar medium per unit mass formed. Yields from isolated stars are
    included only for stars which have evolved off of the main sequence, but the Type Ia supernova contribution is included at
    every mass: it arises from binaries whose secondary has yet to evolve, and is returned per unit interval of secondary mass.
    !!}
    implicit none
    double precision, intent(in   ) :: massInitial

    if (stellarAstrophysics_%lifetime(massInitial,metallicityNode) < ageNode) then
       integrandYield=+initialMassFunction_%phi      (massInitial                ) &
            &         *stellarAstrophysics_%massYield(massInitial,metallicityNode)
    else
       integrandYield=+0.0d0
    end if
    integrandYield   =+integrandYield                                                                       &
         &            +initialMassFunction_%phi  (massInitial                                             ) &
         &            *supernovaeTypeIa_   %yield(initialMassFunction_,massInitial,ageNode,metallicityNode)
    return
  end function integrandYield

end program Test_Stellar_Populations


