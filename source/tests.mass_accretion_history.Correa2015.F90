!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a program which tests the \cite{correa_accretion_2015} halo mass formation history.

program Test_Correa2015_MAH
  !% Tests the \cite{correa_accretion_2015} halo mass formation history algorithm.
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters
  use Dark_Matter_Halo_Mass_Accretion_Histories
  use Cosmology_Functions
  use Unit_Tests
  use Galacticus_Input_Paths
  use Galacticus_Nodes
  implicit none
  type            (treeNode                               ), pointer      :: node
  class           (nodeComponentBasic                     ), pointer      :: basic
  class           (cosmologyFunctionsClass                ), pointer      :: cosmologyFunctions_
  class           (darkMatterHaloMassAccretionHistoryClass), pointer      :: darkMatterHaloMassAccretionHistory_
  type            (varying_string                         )               :: parameterFile
  double precision                                         , dimension(3) :: time                               , mass          , &
       &                                                                     redshift                           , redshiftTarget
  integer                                                                 :: i

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.mass_accretion_history.Correa2015.size')
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Correa et al. 2015 mass accretion history algorithms")
  ! Test Correa et al. 2015 algorithm.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Correa2015.xml'
  call Input_Parameters_File_Open(parameterFile)
  ! Create a node.
  node  => treeNode      (                 )
  ! Get the basic component.
  basic => node    %basic(autoCreate=.true.)
  ! Get required objects.
  cosmologyFunctions_                 => cosmologyFunctions                ()
  darkMatterHaloMassAccretionHistory_ => darkMatterHaloMassAccretionHistory()
  ! Specify halo masses and redshifts.
  mass=          [                  &
       &          1.00000000000d12, &
       &          5.64259219617d11, &
       &          2.96360394164d11  &
       &         ]
  redshiftTarget=[                  &
       &          0.00000000000d00, &
       &          1.00000000000d00, &
       &          2.00000000000d00  &
       &         ]
  ! Set node properties.
  call basic%massSet(                                                                    &
       &                                                               mass          (1) &
       &            )
  call basic%timeSet(                                                                    &
       &             cosmologyFunctions_ %cosmicTime                 (                   &
       &              cosmologyFunctions_%expansionFactorFromRedshift (                  &
       &                                                               redshiftTarget(1) &
       &                                                              )                  &
       &                                                             )                   &
       &            )
  ! Compute the mass accretion history.
  do i=1,size(mass)
     time    (i)=darkMatterHaloMassAccretionHistory_%time(node,mass(i))
     redshift(i)=cosmologyFunctions_ %redshiftFromExpansionFactor(         &
          &       cosmologyFunctions_%expansionFactor             (        &
          &                                                        time(i) &
          &                                                       )        &
          &                                                      )
  end do  
  call Assert('mass accretion history',1.0d0+redshift,1.0d0+redshiftTarget,relTol=1.0d-3)
  ! Close the input parameter file.
  call Input_Parameters_File_Close()
  ! End unit tests.
  call Unit_Tests_End_Group       ()
  call Unit_Tests_Finish          ()

end program Test_Correa2015_MAH
