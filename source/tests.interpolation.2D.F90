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
Contains a program to test 2D interpolation routines.
!!}

program Test_Interpolation_2D
  !!{
  Tests that 2D interpolation routines work.
  !!}
  use :: Display                             , only : displayVerbositySet     , verbosityLevelStandard
  use :: Numerical_Interpolation_2D_Irregular, only : Interpolate_2D_Irregular, interp2dIrregularObject
  use :: Unit_Tests                          , only : Assert                  , Unit_Tests_Begin_Group , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                         , dimension(30) :: xTable                       , yTable    , &
       &                                                      zTable
  double precision                         , dimension(10) :: x                            , y         , &
       &                                                      z                            , zExpected
  type            (interp2dIrregularObject)                :: interpolationWorkspace
  logical                                                  :: resetInterpolation    =.true.
  integer                                                  :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("2D interpolation")

  xTable=[                &
       &  0.5858591716d0, &
       &  0.2593247085d0, &
       &  0.6067972188d0, &
       &  0.9650516119d0, &
       &  0.1565019707d0, &
       &  0.0368830161d0, &
       &  0.7361121764d0, &
       &  0.6132222861d0, &
       &  0.8829957652d0, &
       &  0.1136262016d0, &
       &  0.1526421341d0, &
       &  0.2000633636d0, &
       &  0.5116603314d0, &
       &  0.3130675307d0, &
       &  0.2845182032d0, &
       &  0.1251186002d0, &
       &  0.6374057555d0, &
       &  0.2599673290d0, &
       &  0.3532640981d0, &
       &  0.2434312925d0, &
       &  0.8115300918d0, &
       &  0.7605051459d0, &
       &  0.9223960638d0, &
       &  0.7539049876d0, &
       &  0.6453876849d0, &
       &  0.8891475596d0, &
       &  0.6340141268d0, &
       &  0.9668940594d0, &
       &  0.2184850988d0, &
       &  0.4083093568d0  &
       & ]

  yTable=[                &
       &  0.4029597505d0, &
       &  0.5841232045d0, &
       &  0.9449225385d0, &
       &  0.7491419930d0, &
       &  0.4790344597d0, &
       &  0.8948163614d0, &
       &  0.5134521471d0, &
       &  0.3998105358d0, &
       &  0.4557386842d0, &
       &  0.9514852846d0, &
       &  0.5106885526d0, &
       &  0.1798609016d0, &
       &  0.7271458325d0, &
       &  0.5617648023d0, &
       &  0.8786791167d0, &
       &  0.1111354376d0, &
       &  0.2681871885d0, &
       &  0.5201914045d0, &
       &  0.4110387624d0, &
       &  0.0812157919d0, &
       &  0.6856181254d0, &
       &  0.8174362355d0, &
       &  0.5657368945d0, &
       &  0.1088217734d0, &
       &  0.7209716486d0, &
       &  0.7982971100d0, &
       &  0.8367323740d0, &
       &  0.1535944170d0, &
       &  0.6516215871d0, &
       &  0.3544660839d0  &
       & ]

  zTable=[                &
       &  0.1852220492d0, &
       &  0.1658587228d0, &
       &  0.2764434959d0, &
       &  0.9114234443d0, &
       &  0.1511786538d0, &
       &  0.3403952047d0, &
       &  0.1516256067d0, &
       &  0.5173773007d0, &
       &  0.2840934172d0, &
       &  0.2910639355d0, &
       &  0.3543079421d0, &
       &  0.8754325090d0, &
       &  0.8899622862d0, &
       &  0.7871556818d0, &
       &  0.8994970354d0, &
       &  0.3769161929d0, &
       &  0.1069722278d0, &
       &  0.9743329831d0, &
       &  0.3381532328d0, &
       &  0.7622058052d0, &
       &  0.8541964456d0, &
       &  0.3940045158d0, &
       &  0.4738280624d0, &
       &  0.3373845336d0, &
       &  0.3812024156d0, &
       &  0.3954334166d0, &
       &  0.0755188344d0, &
       &  0.6542130164d0, &
       &  0.9529872043d0, &
       &  0.1231806851d0  &
       & ]

  x=     [                &
       &  0.8616389739d0, &
       &  0.6941142939d0, &
       &  0.6693431344d0, &
       &  0.5921686948d0, &
       &  0.5993416677d0, &
       &  0.8380976762d0, &
       &  0.1836908497d0, &
       &  0.1333316481d0, &
       &  0.8790745125d0, &
       &  0.7046271069d0  &
       & ]

  y=     [                &
       &  0.7180532431d0, &
       &  0.0736140739d0, &
       &  0.5607836703d0, &
       &  0.1566778761d0, &
       &  0.0984902089d0, &
       &  0.2334855772d0, &
       &  0.9065899872d0, &
       &  0.5834639692d0, &
       &  0.3272385504d0, &
       &  0.5767996241d0  &
       & ]

  zExpected=[                       &
       &     0.77761029436913132d0, &
       &     0.26662388093008405d0, &
       &     0.22707976663464668d0, &
       &     0.14645872668613416d0, &
       &     0.27287015778075080d0, &
       &     0.35878910935228742d0, &
       &     0.56940253473405611d0, &
       &     0.41993939199925162d0, &
       &     0.25723639212386112d0, &
       &     0.22254465968375139d0  &
       &    ]

  do i=1,size(x)
     z(i)=Interpolate_2D_Irregular(xTable,yTable,zTable,x(i),y(i),interpolationWorkspace,reset=resetInterpolation)
  end do
  call Assert('random points',z,zExpected,absTol=1.0d-6,relTol=1.0d-6)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Interpolation_2D
