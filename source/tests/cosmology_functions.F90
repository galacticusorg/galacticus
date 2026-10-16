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
Contains a program to test the determinism of the tabulations built by the cosmology functions classes.
!!}

program Test_Cosmology_Functions
  !!{RST
  Tests that the tabulations built by the cosmology functions classes are independent of the sequence of epochs at which they
  happen to be asked for values.

  Each class tabulates expansion factor against cosmic time - and, for the matter plus cosmological constant class, comoving
  distance against cosmic time - extending the tabulation whenever a request falls outside it. With the abscissae pinned to an
  absolute lattice the tabulation reached after a given set of requests depends only on which lattice points those requests
  span, and not on the order in which they arrived. That is asserted here by presenting the same set of epochs to two objects
  in opposite orders and requiring the results to agree *exactly*.

  Every epoch is requested before any value is recorded. This is deliberate: the expansion factor is integrated forward from
  an initial condition imposed at the earliest tabulated epoch, so a tabulation which has not yet been asked for the earliest
  epoch of the set is not the same tabulation, and could not be expected to agree with one which has. What is asserted is that
  the tabulation *arrived at* is the same, which is the property this determinism buys.
  !!}
  use :: Cosmology_Functions , only : cosmologyFunctionsMatterLambda, cosmologyFunctionsMatterDarkEnergy
  use :: Cosmology_Parameters, only : cosmologyParametersSimple
  use :: Display             , only : displayVerbositySet           , verbosityLevelStandard
  use :: Unit_Tests          , only : Assert                        , Unit_Tests_Begin_Group            , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  ! Epochs at which each tabulation is queried, in gigayears. They span several decades, so that every tabulation must be
  ! extended repeatedly, and they are deliberately not points of the lattice. The expansion factor tabulations begin life
  ! spanning 10⁻⁴ to 20 Gyr, so the epochs presented to them must reach beyond that at *both* ends or no extension is provoked
  ! at all - the matter plus dark energy case initially did not, and so passed even against the unpinned tabulation. The
  ! comoving distance tabulation ends at the present day, so the epochs presented to it stop below the present day age of this
  ! cosmology (13.5 Gyr).
  double precision                                    , dimension(7), parameter :: timesExpansion      =[2.1d-5,3.7d-3,4.10d-2,0.83d0,2.9d0,9.1d0,47.0d0]
  double precision                                    , dimension(7), parameter :: timesDistance       =[3.7d-3,4.1d-2,0.83d+0,2.90d0,5.3d0,9.1d0,13.0d0]
  type            (cosmologyParametersSimple         )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda    )                          :: lambdaAscending_    , lambdaDescending_
  type            (cosmologyFunctionsMatterDarkEnergy)                          :: darkEnergyAscending_, darkEnergyDescending_
  double precision                                    , dimension(7)            :: resultAscending     , resultDescending
  integer                                                                       :: i                   , pass

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cosmology functions: tabulation determinism")
  !![
  <referenceConstruct object="cosmologyParameters_" >
   <constructor>
    cosmologyParametersSimple         (                                                  &amp;
     &amp;                             OmegaMatter                = 0.30d0             , &amp;
     &amp;                             OmegaBaryon                = 0.00d0             , &amp;
     &amp;                             OmegaDarkEnergy            = 0.70d0             , &amp;
     &amp;                             temperatureCMB             = 2.78d0             , &amp;
     &amp;                             HubbleConstant             =70.00d0               &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="lambdaAscending_"     >
   <constructor>
    cosmologyFunctionsMatterLambda    (                                                  &amp;
     &amp;                             cosmologyParameters_       =cosmologyParameters_  &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="lambdaDescending_"    >
   <constructor>
    cosmologyFunctionsMatterLambda    (                                                  &amp;
     &amp;                             cosmologyParameters_       =cosmologyParameters_  &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkEnergyAscending_" >
   <constructor>
    cosmologyFunctionsMatterDarkEnergy(                                                  &amp;
     &amp;                             cosmologyParameters_       =cosmologyParameters_, &amp;
     &amp;                             darkEnergyEquationOfStateW0=-0.80d0             , &amp;
     &amp;                             darkEnergyEquationOfStateW1= 0.00d0               &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkEnergyDescending_">
   <constructor>
    cosmologyFunctionsMatterDarkEnergy(                                                  &amp;
     &amp;                             cosmologyParameters_       =cosmologyParameters_, &amp;
     &amp;                             darkEnergyEquationOfStateW0=-0.80d0             , &amp;
     &amp;                             darkEnergyEquationOfStateW1= 0.00d0               &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  !!]

  ! Expansion factor vs. time, matter plus cosmological constant.
  do pass=1,2
     do i=1,size(timesExpansion)
        resultAscending (i)=lambdaAscending_ %expansionFactor(timesExpansion(                       i))
        resultDescending(i)=lambdaDescending_%expansionFactor(timesExpansion(size(timesExpansion)+1-i))
     end do
  end do
  call Assert("matter + cosmological constant: expansion factor is independent of the order in which epochs are requested",resultAscending,reverse(resultDescending),absTol=0.0d0)

  ! Cosmic time vs. expansion factor - the inverse of the same tabulation.
  do pass=1,2
     do i=1,size(timesDistance)
        resultAscending (i)=lambdaAscending_ %cosmicTime(lambdaAscending_ %expansionFactor(timesDistance(                      i)))
        resultDescending(i)=lambdaDescending_%cosmicTime(lambdaDescending_%expansionFactor(timesDistance(size(timesDistance)+1-i)))
     end do
  end do
  call Assert("matter + cosmological constant: cosmic time is independent of the order in which epochs are requested"     ,resultAscending,reverse(resultDescending),absTol=0.0d0)

  ! Comoving distance vs. time.
  do pass=1,2
     do i=1,size(timesDistance)
        resultAscending (i)=lambdaAscending_ %distanceComoving(timesDistance(                      i))
        resultDescending(i)=lambdaDescending_%distanceComoving(timesDistance(size(timesDistance)+1-i))
     end do
  end do
  call Assert("matter + cosmological constant: comoving distance is independent of the order in which epochs are requested",resultAscending,reverse(resultDescending),absTol=0.0d0)

  ! Expansion factor vs. time, matter plus dark energy.
  do pass=1,2
     do i=1,size(timesExpansion)
        resultAscending (i)=darkEnergyAscending_ %expansionFactor(timesExpansion(                       i))
        resultDescending(i)=darkEnergyDescending_%expansionFactor(timesExpansion(size(timesExpansion)+1-i))
     end do
  end do
  call Assert("matter + dark energy: expansion factor is independent of the order in which epochs are requested"           ,resultAscending,reverse(resultDescending),absTol=0.0d0)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  function reverse(values)
    !!{RST
    Return ``values`` in reverse order.
    !!}
    implicit none
    double precision, intent(in   ), dimension(           :) :: values
    double precision               , dimension(size(values)) :: reverse
    integer                                                  :: i

    do i=1,size(values)
       reverse(i)=values(size(values)-i+1)
    end do
    return
  end function reverse

end program Test_Cosmology_Functions
