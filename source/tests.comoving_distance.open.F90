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

program Tests_Comoving_Distance_Open
  !% Tests comoving distance calculations for an open Universe. Distances calculated using Python implementation of Ned Wright's
  !% cosmology calculator available from: http://www.astro.ucla.edu/~wright/CC.python
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Cosmology_Functions
  use Cosmology_Functions_Options
  use Memory_Management
  implicit none
  double precision                , dimension(8), parameter :: redshift     =[0.1d0,1.0d0,3.0d0,9.0d0,30.0d0,100.0d0,300.0d0,1000.0d0]                                                            
  double precision                , dimension(8), parameter :: distanceOpen =[2.836546d0,19.566092d0,36.248716d0,53.313669d0,67.286209d0,75.856732d0,80.401477d0,83.217159d0]                     
  type            (varying_string)                          :: parameterFile                                                                                                                      
  character       (len=1024      )                          :: message                                                                                                                            
  integer                                                   :: iExpansion                                                                                                                         
  double precision                                          :: distance                                                                                                      , distanceModulus, & 
       &                                                       time                                                                                                          , timeLookup         
  
  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.comoving_distance.open.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Comoving distance: open cosmology")

  ! Test comoving distance in an open universe.
  parameterFile='testSuite/parameters/comovingDistance/open.xml'
  call Input_Parameters_File_Open(parameterFile)
  do iExpansion=1,size(redshift)
     time=Cosmology_Age(Expansion_Factor_From_Redshift(redshift(iExpansion)))
     distance=Comoving_Distance(time)
     write (message,'(a,f6.1,a)') "comoving distance [z=",redshift(iExpansion),"]"
     call Assert(trim(message),distance,distanceOpen(iExpansion),relTol=1.0d-3)
     timeLookup=Time_From_Comoving_Distance(distance)
     write (message,'(a,f6.1,a)') "cosmic time [distance=",distance," Mpc]"
     call Assert(trim(message),timeLookup,time,relTol=1.0d-3)
     distance=Comoving_Distance_Conversion(distanceTypeComoving,redshift=redshift(iExpansion))
     write (message,'(a,f6.1,a)') "comoving distance [direct; z=",redshift(iExpansion),"]"
     call Assert(trim(message),distance,distanceOpen(iExpansion),relTol=1.0d-3)
     distanceModulus=25.0d0+5.0d0*log10(distance*(1.0d0+redshift(iExpansion)))
     distance=Comoving_Distance_Conversion(distanceTypeComoving,distanceModulus=distanceModulus)
     write (message,'(a,f6.1,a)') "comoving distance [D=",distanceModulus,"]"
     call Assert(trim(message),distance,distanceOpen(iExpansion),relTol=1.0d-3)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Comoving_Distance_Open
