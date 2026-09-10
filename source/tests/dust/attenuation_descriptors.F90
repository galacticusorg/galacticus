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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program to test the dust attenuation descriptor data structures.
!!}

program Test_Dust_Attenuation_Descriptors
  !!{RST
  Tests the data structures which mediate between luminosity-producing objects and dust attenuation objects: the age
  binning negotiated through a ``decompositionRequest``, and the reduction of an attenuated
  ``luminosityDecomposition`` back into output element values.
  !!}
  use :: Display                     , only : displayVerbositySet , verbosityLevelStandard
  use :: Dust_Attenuation_Descriptors, only : decompositionRequest, emissionDescriptor    , luminosityDecomposition
  use :: Galactic_Structure_Options  , only : componentTypeUnknown
  use :: Unit_Tests                  , only : Assert              , Unit_Tests_Begin_Group, Unit_Tests_End_Group   , Unit_Tests_Finish
  implicit none
  type            (decompositionRequest   )                            :: request
  type            (luminosityDecomposition)                            :: decomposition
  type            (emissionDescriptor     )                            :: descriptor
  double precision                         , allocatable, dimension(:) :: values
  double precision                                                     :: ageMinimum   , ageMaximum

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Dust attenuation descriptors")

  ! An emission descriptor which has not been filled in should describe emission of unknown origin spanning all ages.
  call Unit_Tests_Begin_Group("default emission descriptor")
  call Assert("wavelength unresolved" ,descriptor%wavelength    < 0.0d0                ,.true.)
  call Assert("metallicity unresolved",descriptor%metallicity   < 0.0d0                ,.true.)
  call Assert("radius unresolved"     ,descriptor%radius        < 0.0d0                ,.true.)
  call Assert("age spans all"         ,descriptor%ageMinimum                           ,0.0d0 )
  call Assert("age unbounded above"   ,descriptor%ageMaximum    == huge(0.0d0)         ,.true.)
  call Assert("component unknown"     ,descriptor%componentType == componentTypeUnknown,.true.)
  call Unit_Tests_End_Group()

  ! A request with no age boundaries requests no age resolution at all.
  call Unit_Tests_Begin_Group("age binning: unresolved")
  call Assert("single bin"               ,request%countAgeBins(       ),1     )
  call Assert("all ages in bin 1 (young)",request%ageBinIndex (0.001d0),1     )
  call Assert("all ages in bin 1 (old)"  ,request%ageBinIndex (1.000d1),1     )
  call request%ageBinRange(1,ageMinimum,ageMaximum)
  call Assert("bin starts at zero"       ,ageMinimum                   ,0.0d0 )
  call Assert("bin unbounded above"      ,ageMaximum == huge(0.0d0)    ,.true.)
  call Unit_Tests_End_Group()

  ! A request with two internal boundaries defines three bins. Bins are closed below and open above, so an age
  ! exactly on a boundary falls into the older bin.
  call Unit_Tests_Begin_Group("age binning: two boundaries")
  request%ageBoundaries=[1.0d-2,1.0d-1]
  call Assert("three bins"             ,request%countAgeBins(      ),3)
  call Assert("below first boundary"   ,request%ageBinIndex (5.0d-3),1)
  call Assert("on first boundary"      ,request%ageBinIndex (1.0d-2),2)
  call Assert("between boundaries"     ,request%ageBinIndex (5.0d-2),2)
  call Assert("on second boundary"     ,request%ageBinIndex (1.0d-1),3)
  call Assert("above second boundary"  ,request%ageBinIndex (1.0d+1),3)
  call request%ageBinRange(2,ageMinimum,ageMaximum)
  call Assert("middle bin lower bound" ,ageMinimum,1.0d-2)
  call Assert("middle bin upper bound" ,ageMaximum,1.0d-1)
  call request%ageBinRange(3,ageMinimum,ageMaximum)
  call Assert("upper bin lower bound"  ,ageMinimum,1.0d-1)
  call Assert("upper bin unbounded"    ,ageMaximum == huge(0.0d0),.true.)
  call Unit_Tests_End_Group()

  ! Reduction of a decomposition sums each parcel, weighted by its transmission, into the output element which it
  ! contributes to.
  call Unit_Tests_Begin_Group("reduction")
  call decomposition%initialize(4,2)
  call Assert("term count",decomposition%countTerms(),4)
  decomposition%luminosities=[1.0d0,2.0d0,3.0d0,4.0d0]
  decomposition%elementIndex=[1    ,1    ,2    ,2    ]
  ! With unit transmission the reduction simply sums the parcels of each element.
  call decomposition%reduce([1.0d0,1.0d0,1.0d0,1.0d0],values)
  call Assert("unattenuated element 1",values(1),3.0d+0,relTol=1.0d-12)
  call Assert("unattenuated element 2",values(2),7.0d+0,relTol=1.0d-12)
  ! With a transmission which differs between parcels of the same element, the parcels are weighted individually.
  call decomposition%reduce([1.0d0,5.0d-1,2.5d-1,0.0d0],values)
  call Assert("attenuated element 1"  ,values(1),2.0d+0,relTol=1.0d-12)
  call Assert("attenuated element 2"  ,values(2),7.5d-1,relTol=1.0d-12)
  ! Zero transmission extinguishes everything.
  call decomposition%reduce([0.0d0,0.0d0,0.0d0,0.0d0],values)
  call Assert("fully extinguished",values,[0.0d0,0.0d0])
  call Unit_Tests_End_Group()

  ! An empty decomposition is legal, and reduces to no elements.
  call Unit_Tests_Begin_Group("empty decomposition")
  call decomposition%initialize(0,0)
  call Assert("no terms",decomposition%countTerms(),0)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Dust_Attenuation_Descriptors
