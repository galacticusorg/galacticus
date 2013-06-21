!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
!% \href{http://cosmo.nyu.edu/~tinker/massfunction/MF_code.tar}{code}.

program Tests_Halo_Mass_Function_Tinker
  !% Tests the \cite{tinker_towardhalo_2008} mass function by comparing to Jeremy Tinker's
  !% \href{http://cosmo.nyu.edu/~tinker/massfunction/MF_code.tar}{code}.
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Memory_Management
  use Halo_Mass_Function
  use Cosmology_Functions
  use Cosmological_Parameters
  use File_Utilities
  use Critical_Overdensity
  implicit none
  type            (varying_string)                            :: parameterFile                                    
  integer                                                     :: fUnit        , i           , massCount           
  double precision                                            :: time                                             
  double precision                , allocatable, dimension(:) :: mass         , massFunction, massFunctionTinker  
  logical                         , allocatable, dimension(:) :: success                                          
  
  ! Read in basic code memory usage.                                                                                                             
  call Code_Memory_Usage('tests.halo_mass_function.Tinker.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Halo mass function: Tinker et al. (2008)")

  ! Test the Tinker et al. (2008) dark matter halo mass function.
  parameterFile='testSuite/parameters/haloMassFunction/tinker.xml'
  call Input_Parameters_File_Open(parameterFile)
  time=Cosmology_Age(1.0d0)

  ! Determine number of masses in reference data file and allocate arrays.
  massCount=Count_Lines_In_File('testSuite/data/haloMassFunction/tinker.txt')
  call Alloc_Array(mass              ,[massCount])
  call Alloc_Array(massFunction      ,[massCount])
  call Alloc_Array(massFunctionTinker,[massCount])
  call Alloc_Array(success           ,[massCount])

  ! Ensure that critical density and critical overdensity for collapse are consistent with values used in our input file to
  ! Tinker's code.
  call Assert('critical density consistency'                 ,Critical_Density()/Little_H_0()**2     ,2.7751950000000000d11,relTol=1.0d-6)
  call Assert('critical overdensity for collapse consistency',Critical_Overdensity_for_Collapse(time),1.6755779626281502d00,relTol=1.0d-6)

  ! Compute mass function for each reference mass.
  open(newUnit=fUnit,file='testSuite/data/haloMassFunction/tinker.txt',status='old',form='formatted')
  do i=1,massCount
     read (fUnit,*) mass(i),massFunctionTinker(i)
     mass              (i)=mass              (i)/Little_H_0()
     massFunctionTinker(i)=massFunctionTinker(i)*Little_H_0()**4
     massFunction      (i)=Halo_Mass_Function_Differential(time,mass(i))
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
  call Input_Parameters_File_Close

  ! Clean up memory.
  call Dealloc_Array(mass              )
  call Dealloc_Array(massFunction      )
  call Dealloc_Array(massFunctionTinker)
  call Dealloc_Array(success           )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Tests_Halo_Mass_Function_Tinker
