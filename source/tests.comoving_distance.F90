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

program Tests_Comoving_Distance
  !% Tests comoving distance calculations for various universes. Distances calculated using Python
  !% \href{http://www.astro.ucla.edu/~wright/CC.python}{implementation} of Ned Wright's cosmology
  !% calculator.
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Cosmology_Functions_Options
  use Cosmology_Functions
  use Cosmology_Parameters
  use Memory_Management
  implicit none
  double precision                                , dimension(8), parameter         :: redshift                               =[0.1000000d0,1.0000000d0,3.0000000d0,9.0000000d0,30.0000000d0,100.0000000d0,300.0000000d0,1000.0000000d0]
  double precision                                , dimension(8)           , target :: distanceEdS                            =[27.9031290d0,175.6143280d0,299.7923450d0,409.9791330d0,491.8947530d0,539.9166970d0,564.9913370d0,580.4490760d0]
  double precision                                , dimension(8)           , target :: distanceOpen                           =[2.8365460d0,19.5660920d0,36.2487160d0,53.3136690d0,67.2862090d0,75.8567320d0,80.4014770d0,83.2171590d0]
  double precision                                , dimension(8)           , target :: distanceCosmologicalConstant           =[2.9291813d0,23.1267935d0,44.4897529d0,64.4722404d0,79.4221370d0,88.1893579d0,92.7669523d0,95.5884181d0]
  double precision                                , dimension(:), pointer           :: thisDistance
  type            (cosmologyParametersSimple     )                                  :: cosmologyParametersCosmologicalConstant                                                                                                                 , cosmologyParametersOpen
  class           (cosmologyFunctionsClass       )              , pointer           :: cosmologyFunctionsEdS                                                                                                                                   , thisCosmologyFunctions
  type            (cosmologyFunctionsMatterLambda)                         , target :: cosmologyFunctionsCosmologicalConstant                                                                                                                  , cosmologyFunctionsOpen
  type            (varying_string                )                                  :: parameterFile
  character       (len=1024                      )                                  :: message
  integer                                                                           :: i                                                                                                                                                       , iExpansion
  double precision                                                                  :: distance                                                                                                                                                , distanceModulus        , &
       &                                                                               time                                                                                                                                                    , timeLookup

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.comoving_distance.size')
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Comoving distance")
  ! Cosmology functions for in an Einstein-de Sitter universe. For this case, we use the default settings.
  parameterFile='testSuite/parameters/comovingDistance/EdS.xml'
  call Input_Parameters_File_Open(parameterFile)
  cosmologyFunctionsEdS => cosmologyFunctions()
  ! Cosmology functions for an open Universe.
  cosmologyParametersOpen                =cosmologyParametersSimple     (0.3d0,0.0d0,0.0d0,2.78d0,10000.0d0     )
  cosmologyFunctionsOpen                 =cosmologyFunctionsMatterLambda(cosmologyParametersOpen                )
  ! Cosmology functions for a cosmological constant Universe.
  cosmologyParametersCosmologicalConstant=cosmologyParametersSimple     (0.3d0,0.0d0,0.7d0,2.78d0,10000.0d0     )
  cosmologyFunctionsCosmologicalConstant =cosmologyFunctionsMatterLambda(cosmologyParametersCosmologicalConstant)
  ! Iterate over cosmologies.
  do i=1,3
     select case (i)
     case (1)
        call Unit_Tests_Begin_Group("Einstein-de Sitter"   )
        thisCosmologyFunctions => cosmologyFunctionsEdS
        thisDistance           => distanceEdS
     case (2)
        call Unit_Tests_Begin_Group("Open universe"        )
        thisCosmologyFunctions => cosmologyFunctionsOpen
        thisDistance           => distanceOpen
     case (3)
        call Unit_Tests_Begin_Group("Cosmological constant")
        thisCosmologyFunctions => cosmologyFunctionsCosmologicalConstant
        thisDistance           => distanceCosmologicalConstant
     end select
     do iExpansion=1,size(redshift)
        time=thisCosmologyFunctions%cosmicTime(thisCosmologyFunctions%expansionFactorFromRedshift(redshift(iExpansion)))
        distance=thisCosmologyFunctions%distanceComoving(time)
        write (message,'(a,f6.1,a)') "comoving distance [        z=",redshift(iExpansion),"    ]"
        call Assert(trim(message),distance,thisDistance(iExpansion),relTol=1.0d-3)
        timeLookup=thisCosmologyFunctions%timeAtDistanceComoving(distance)
        write (message,'(a,f6.1,a)') "cosmic time       [distance =",distance            ," Mpc]"
        call Assert(trim(message),timeLookup,time,relTol=1.0d-3)
        distance=thisCosmologyFunctions%distanceComovingConvert(distanceTypeComoving,redshift=redshift(iExpansion))
        write (message,'(a,f6.1,a)') "comoving distance [direct; z=",redshift(iExpansion),"    ]"
        call Assert(trim(message),distance,thisDistance(iExpansion),relTol=1.0d-3)
        distanceModulus=25.0d0+5.0d0*log10(distance*(1.0d0+redshift(iExpansion)))
        distance=thisCosmologyFunctions%distanceComovingConvert(distanceTypeComoving,distanceModulus=distanceModulus)
        write (message,'(a,f6.1,a)') "comoving distance [        D=",distanceModulus     ,"    ]"
        call Assert(trim(message),distance,thisDistance(iExpansion),relTol=1.0d-3)
     end do
     call Unit_Tests_End_Group()
  end do
  call Input_Parameters_File_Close()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Comoving_Distance
