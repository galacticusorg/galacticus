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

program Tests_Linear_Growth_EdS
  !!{RST
  Tests linear growth calculations.
  !!}
  use :: Cosmology_Parameters, only : cosmologyParametersSimple
  use :: Cosmology_Functions , only : cosmologyFunctionsMatterLambda
  use :: Display             , only : displayVerbositySet           , verbosityLevelStandard
  use :: Linear_Growth       , only : componentDarkMatter           , linearGrowthCollisionlessMatter
  use :: Unit_Tests          , only : Assert                        , Unit_Tests_Begin_Group         , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                                 , dimension(8), parameter :: redshift                =[0.0d0,1.0d0,3.0d0,9.0d0,30.0d0,100.0d0,300.0d0,1000.0d0]
  type            (cosmologyParametersSimple      )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda )                          :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter)                          :: linearGrowth_
  character       (len=1024                       )                          :: message
  integer                                                                    :: iExpansion
  double precision                                                           :: expansionFactor                                                                   , linearGrowthFactor
  ! An epoch beyond any reached above, used to force the tabulation to be extended.
  double precision                                               , parameter :: expansionFactorExtension=4.0d0
  double precision                                 , dimension(8)            :: linearGrowthFactorBefore                                                          , linearGrowthFactorAfter
  logical                                                                    :: unchangedByExtension

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Linear growth: Einstein-de Sitter")
  ! Test growth factor in an Einstein-de Sitter universe. Growth factor should equal the expansion factor.
  !![
  <referenceConstruct object="cosmologyParameters_">
   <constructor>
    cosmologyParametersSimple      (                                           &amp;
     &amp;                          OmegaMatter         = 1.00d0             , &amp;
     &amp;                          OmegaBaryon         = 0.00d0             , &amp;
     &amp;                          OmegaDarkEnergy     = 0.00d0             , &amp;
     &amp;                          temperatureCMB      = 2.78d0             , &amp;
     &amp;                          HubbleConstant      =73.00d0               &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_" >
   <constructor>
    cosmologyFunctionsMatterLambda (                                           &amp;
     &amp;                          cosmologyParameters_=cosmologyParameters_  &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"       >
   <constructor>
    linearGrowthCollisionlessMatter(                                           &amp;
     &amp;                          cosmologyParameters_=cosmologyParameters_, &amp;
     &amp;                          cosmologyFunctions_ =cosmologyFunctions_   &amp;
     &amp;                         )
   </constructor>
  </referenceConstruct>
  !!]
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(redshift(iExpansion))
     linearGrowthFactor=linearGrowth_%value(expansionFactor=expansionFactor,component=componentDarkMatter)
     write (message,'(a,f6.1,a)') "dark matter linear growth factor [z=",redshift(iExpansion),"]"
     call Assert(trim(message),linearGrowthFactor,expansionFactor,relTol=1.0d-3)
  end do
  ! Check that extending the tabulation leaves the values it already held unchanged. The tabulated epochs lie on an absolute
  ! lattice, so a request beyond the latest of them adds points to the tabulation and resumes the integration of the growth
  ! factor from the state in which it was left, rather than laying a fresh grid across the enlarged range and integrating it
  ! again from the start - so the growth factor at every epoch already tabulated must be reproduced bit for bit. The values are
  ! recorded only after the loop above, so that the comparison is across the extension alone: the tabulation reaches its final
  ! range during that loop, and a request for an epoch *earlier* than any tabulated does rebuild it - the growth factor is
  ! integrated forward from the earliest tabulated epoch, so moving that epoch necessarily changes every value.
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(redshift(iExpansion))
     linearGrowthFactorBefore(iExpansion)=linearGrowth_%value(expansionFactor=expansionFactor,component=componentDarkMatter)
  end do
  linearGrowthFactor=linearGrowth_%value(expansionFactor=expansionFactorExtension,component=componentDarkMatter)
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(redshift(iExpansion))
     linearGrowthFactorAfter(iExpansion)=linearGrowth_%value(expansionFactor=expansionFactor,component=componentDarkMatter)
  end do
  unchangedByExtension=all(linearGrowthFactorAfter == linearGrowthFactorBefore)
  call Assert("dark matter linear growth factor is unchanged by extension of the tabulation",unchangedByExtension,.true.)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Linear_Growth_EdS
