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

!% Contains a program which tests the \cite{tinker_towardhalo_2008} mass function by comparing to Jeremy Tinker's
!% \href{https://cosmo.nyu.edu/~tinker/massfunction/MF_code.tar}{code}.

program Tests_Halo_Mass_Function_Tinker
  !% Tests the \cite{tinker_towardhalo_2008} mass function by comparing to Jeremy Tinker's
  !% \href{https://cosmo.nyu.edu/~tinker/massfunction/MF_code.tar}{code}.
  use :: Cosmological_Density_Field, only : criticalOverdensity           , criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctions            , cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParameters           , cosmologyParametersClass, hubbleUnitsLittleH
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: File_Utilities            , only : Count_Lines_In_File
  use :: Galacticus_Display        , only : Galacticus_Verbosity_Level_Set, verbosityStandard
  use :: Halo_Mass_Functions       , only : haloMassFunction              , haloMassFunctionClass
  use :: ISO_Varying_String        , only : varying_string                , assignment(=)
  use :: Input_Parameters          , only : inputParameters
  use :: Memory_Management         , only : allocateArray                 , deallocateArray
  use :: Unit_Tests                , only : Assert                        , Unit_Tests_Begin_Group  , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (varying_string          )                                     :: parameterFile
  integer                                                                        :: fUnit                    , i           , &
       &                                                                            massCount
  double precision                                                               :: time
  double precision                                   , allocatable, dimension(:) :: mass                     , massFunction, &
       &                                                                            massFunctionTinker
  logical                                            , allocatable, dimension(:) :: success
  class           (cosmologyParametersClass), pointer                            :: cosmologyParameters_
  class           (cosmologyFunctionsClass ), pointer                            :: cosmologyFunctions_
  class           (criticalOverdensityClass), pointer                            :: criticalOverdensity_
  class           (haloMassFunctionClass   ), pointer                            :: haloMassFunction_
  type            (inputParameters         )                                     :: parameters

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Halo mass function: Tinker et al. (2008)")

  ! Test the Tinker et al. (2008) dark matter halo mass function.
  parameterFile='testSuite/parameters/haloMassFunction/tinker.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Get required objects.
  cosmologyParameters_ => cosmologyParameters()
  cosmologyFunctions_  => cosmologyFunctions ()
  criticalOverdensity_ => criticalOverdensity()
  haloMassFunction_    => haloMassFunction   ()

  time=cosmologyFunctions_%cosmicTime(1.0d0)

  ! Determine number of masses in reference data file and allocate arrays.
  massCount=Count_Lines_In_File('testSuite/data/haloMassFunction/tinker.txt')
  call allocateArray(mass              ,[massCount])
  call allocateArray(massFunction      ,[massCount])
  call allocateArray(massFunctionTinker,[massCount])
  call allocateArray(success           ,[massCount])

  ! Ensure that critical density and critical overdensity for collapse are consistent with values used in our input file to
  ! Tinker's code.
  call Assert('critical density consistency'                 ,cosmologyParameters_%densityCritical(    )/cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2,2.7751950000000000d11,relTol=1.0d-6)
  call Assert('critical overdensity for collapse consistency',criticalOverdensity_%value          (time)                                                           ,1.6755779626281502d00,relTol=1.0d-6)

  ! Compute mass function for each reference mass.
  open(newUnit=fUnit,file='testSuite/data/haloMassFunction/tinker.txt',status='old',form='formatted')
  do i=1,massCount
     read (fUnit,*) mass(i),massFunctionTinker(i)
     mass              (i)=mass              (i)/cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
     massFunctionTinker(i)=massFunctionTinker(i)*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**4
     massFunction      (i)=haloMassFunction_%differential(time,mass(i))
  end do

  ! Assert that our mass function agrees with the reference data.
  success=                                                       &
       &   abs(1.0d0-massFunction/massFunctionTinker) < 1.0d-02  &
       &  .or.                                                   &
       &   abs(      massFunction-massFunctionTinker) < 1.0d-20
  call Assert(                           &
       &      'halo mass mass function', &
       &      all(success),              &
       &      .true.                     &
       &     )

  ! Clean up memory.
  call deallocateArray(mass              )
  call deallocateArray(massFunction      )
  call deallocateArray(massFunctionTinker)
  call deallocateArray(success           )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Tests_Halo_Mass_Function_Tinker
