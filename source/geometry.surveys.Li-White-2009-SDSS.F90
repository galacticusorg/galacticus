!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the survey geometry used by \cite{li_distribution_2009}.

module Geometry_Surveys_Li_White_2009_SDSS
  !% Implements the survey geometry used by \cite{li_distribution_2009}.
  implicit none
  private
  public :: Geometry_Surveys_Li_White_2009_SDSS_Initialize

contains

  !# <surveyGeometryMethod>
  !#  <unitName>Geometry_Surveys_Li_White_2009_SDSS_Initialize</unitName>
  !# </surveyGeometryMethod>
  subroutine Geometry_Surveys_Li_White_2009_SDSS_Initialize(surveyGeometryMethod,Geometry_Survey_Distance_Maximum_Get&
       &,Geometry_Survey_Solid_Angle_Get ,Geometry_Survey_Volume_Maximum_Get,Geometry_Survey_Window_Functions_Get)
    !% Initializes the ``Li-White-2009-SDSS'' survey geometry module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string  ),          intent(in   ) :: surveyGeometryMethod
    procedure(double precision), pointer, intent(inout) :: Geometry_Survey_Distance_Maximum_Get,Geometry_Survey_Solid_Angle_Get&
         &,Geometry_Survey_Volume_Maximum_Get
    procedure(                ), pointer, intent(inout) :: Geometry_Survey_Window_Functions_Get

    if (surveyGeometryMethod == 'Li-White-2009-SDSS') then
       Geometry_Survey_Distance_Maximum_Get => Geometry_Survey_Distance_Maximum_Li_White_2009_SDSS
       Geometry_Survey_Solid_Angle_Get      => Geometry_Survey_Solid_Angle_Li_White_2009_SDSS
       Geometry_Survey_Volume_Maximum_Get   => Geometry_Survey_Volume_Maximum_Li_White_2009_SDSS
       Geometry_Survey_Window_Functions_Get => Geometry_Survey_Window_Functions_Li_White_2009_SDSS
    end if
    return
  end subroutine Geometry_Surveys_Li_White_2009_SDSS_Initialize

  double precision function Geometry_Survey_Distance_Maximum_Li_White_2009_SDSS(mass)
    !% Compute the maximum distance at which a galaxy is visible.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    implicit none
    double precision, intent(in) :: mass
    double precision             :: redshift,logarithmicMass
    
    ! Find the limiting redshift for this mass using a fit derived from Millennium Simulation SAMs. (See
    ! constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07/massLuminosityRelation.pl for details.)
    logarithmicMass=log10(mass)
    redshift=-5.9502006195004d0+logarithmicMass*(2.63793788603951d0+logarithmicMass*(-0.421075858899237d0+logarithmicMass&
         &*(0.0285198776926787d0+logarithmicMass*(-0.000678327494720407d0))))
    
    ! Convert from redshift to comoving distance.
    Geometry_Survey_Distance_Maximum_Li_White_2009_SDSS=Comoving_Distance_Conversion(output=distanceTypeComoving,redshift&
         &=redshift)
    return
  end function Geometry_Survey_Distance_Maximum_Li_White_2009_SDSS

  double precision function Geometry_Survey_Solid_Angle_Li_White_2009_SDSS()
    !% Return the solid angle of the \cite{li_distribution_2009} sample.
    implicit none
    double precision, parameter :: solidAngleSurvey=2.1901993d0 ! From Percival et al. (2010; MNRAS; 401; 2148)
    
    Geometry_Survey_Solid_Angle_Li_White_2009_SDSS=solidAngleSurvey
    return
  end function Geometry_Survey_Solid_Angle_Li_White_2009_SDSS

  double precision function Geometry_Survey_Volume_Maximum_Li_White_2009_SDSS(mass)
    !% Compute the maximum volume in which a galaxy of given mass could have been observed.
    implicit none
    double precision, intent(in) :: mass

    ! Find the volume associated with this maximum distance.
    Geometry_Survey_Volume_Maximum_Li_White_2009_SDSS=Geometry_Survey_Solid_Angle_Li_White_2009_SDSS()&
         &*Geometry_Survey_Distance_Maximum_Li_White_2009_SDSS(mass)**3/3.0d0    
    return
  end function Geometry_Survey_Volume_Maximum_Li_White_2009_SDSS

  subroutine Geometry_Survey_Window_Functions_Li_White_2009_SDSS(mass1,mass2,boxLength,gridCount,windowFunction1,windowFunction2)
    !% Compute the window function for the survey.
    use FFTW3
    use Vectors
    use Pseudo_Random
    use File_Utilities
    use FGSL
    use Meshes
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Memory_Management
    use Galacticus_Display
    use ISO_Varying_String
    use Cosmology_Functions
    use String_Handling
    use File_Utilities
    use Galacticus_Input_Paths
    use System_Command
    use Galacticus_Error
    implicit none
    double precision         , intent(in   )                                               :: mass1,mass2
    integer                  , intent(in   )                                               :: gridCount
    double precision         , intent(  out)                                               :: boxLength
    complex(c_double_complex), intent(  out),     dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2
    double precision         ,                    dimension(3                            ) :: origin,position1,position2
    integer                                                                                :: i,j
    double precision                                                                       :: comovingDistanceMaximum1&
         &,comovingDistanceMaximum2,comovingDistanceMinimum1,comovingDistanceMinimum2
 
    type   (c_ptr           )                                                              :: plan
    complex(c_double_complex),                    dimension(gridCount,gridCount,gridCount) :: selectionFunction1,selectionFunction2
    complex(c_double_complex)                                                              :: normalization
    logical                  , save                                                        :: geometryInitialized=.false.
    double precision                                                                       :: rightAscension,declination,distance1,distance2
    double precision         , save, allocatable, dimension(:                            ) :: randomTheta,randomPhi
    integer                  , save                                                        :: randomsCount
    integer                                                                                :: randomUnit
    type   (fgsl_rng        ), save                                                        :: pseudoSequenceObject
    logical                  , save                                                        :: reset=.true.
    type   (varying_string  )                                                              :: message
    double precision         , save                                                        :: surveyDistanceMinimum&
         &,surveyDistanceMaximum

    ! Initialize geometry if necessary.
    if (.not.geometryInitialized) then
       ! Randoms file obtained from:  http://sdss.physics.nyu.edu/lss/dr72/random/
       if (.not.File_Exists(Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")) then
          call System_Command_Do("mkdir -p "//Galacticus_Input_Path()//"data/surveyGeometry")
          call System_Command_Do("wget http://sdss.physics.nyu.edu/lss/dr72/random/lss_random-0.dr72.dat -O "//Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")
          if (.not.File_Exists(Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")) call Galacticus_Error_Report('Geometry_Survey_Window_Functions_Li_White_2009_SDSS','unable to download SDSS survey geometry randoms file')
       end if
       randomsCount=Count_Lines_In_File("lss_random-0.dr72.dat")
       call Alloc_Array(randomTheta,[randomsCount])
       call Alloc_Array(randomPhi  ,[randomsCount])
       open(newUnit=randomUnit,file="lss_random-0.dr72.dat",status="old",form="formatted")
       j=0
       do i=1,randomsCount
          read (randomUnit,*) rightAscension,declination
          if     (                                &
               &          rightAscension  > 100.0 &
               &  .and.   rightAscension  < 300.0 &
               &  .and. (                         &
               &          rightAscension  < 247.0 &
               &         .or. declination <  51.0 &
               &        )                         &
               & ) then
             j=j+1
             randomTheta(j)=Pi*(90.0d0-declination   )/180.0d0
             randomPhi  (j)=Pi*        rightAscension /180.0d0
          end if
       end do
       close(randomUnit)
       message="Read "
       message=message//randomsCount//" random points and kept "//j//" of them"
       call Galacticus_Display_Message(message)
       randomsCount=j

       ! Compute the distances corresponding to the minimum and maximum redshifts.
       surveyDistanceMinimum=Comoving_Distance(Cosmology_Age(Expansion_Factor_From_Redshift(0.001d0)))
       surveyDistanceMaximum=Comoving_Distance(Cosmology_Age(Expansion_Factor_From_Redshift(0.500d0)))
       geometryInitialized=.true.
    end if

    ! Find the comoving distance corresponding to this distance module.
    comovingDistanceMaximum1=min(Geometry_Survey_Distance_Maximum_Li_White_2009_SDSS(mass1),surveyDistanceMaximum)
    comovingDistanceMaximum2=min(Geometry_Survey_Distance_Maximum_Li_White_2009_SDSS(mass2),surveyDistanceMaximum)
    comovingDistanceMinimum1=surveyDistanceMinimum
    comovingDistanceMinimum2=surveyDistanceMinimum

    ! Determine a suitable box length for the window function calculation.
    boxLength=3.0d0*max(comovingDistanceMaximum1,comovingDistanceMaximum2)

    ! Set up origin.
    origin=([0.0d0,0.5d0,0.5d0]/dble(gridCount)+[0.5d0,0.5d0,0.5d0])*dble(boxLength)

    ! Populate the cube with the survey selection function.
    selectionFunction1=cmplx(0.0d0,0.0d0)
    selectionFunction2=cmplx(0.0d0,0.0d0)

    ! Loop over randoms.
    do i=1,randomsCount
       ! Choose random distances.
       distance1=(Pseudo_Random_Get(pseudoSequenceObject,reset)*(comovingDistanceMaximum1**3-comovingDistanceMinimum1**3)+comovingDistanceMinimum1**3)**(1.0d0/3.0d0)
       distance2=(Pseudo_Random_Get(pseudoSequenceObject,reset)*(comovingDistanceMaximum2**3-comovingDistanceMinimum2**3)+comovingDistanceMinimum2**3)**(1.0d0/3.0d0)
       ! Convert to Cartesian coordinates.
       position1=distance1*[sin(randomTheta(i))*cos(randomPhi(i)),sin(randomTheta(i))*sin(randomPhi(i)),cos(randomTheta(i))]+origin
       position2=distance2*[sin(randomTheta(i))*cos(randomPhi(i)),sin(randomTheta(i))*sin(randomPhi(i)),cos(randomTheta(i))]+origin
       ! Apply to the mesh.
       call Meshes_Apply_Point(selectionFunction1,boxLength,position1,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
       call Meshes_Apply_Point(selectionFunction2,boxLength,position2,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
    end do
    ! Take the Fourier transform of the selection function.
    plan=fftw_plan_dft_3d(gridCount,gridCount,gridCount,selectionFunction1,windowFunction1,FFTW_FORWARD,FFTW_ESTIMATE)
    call fftw_execute_dft(plan,selectionFunction1,windowFunction1)
    call fftw_execute_dft(plan,selectionFunction2,windowFunction2)
    call fftw_destroy_plan(plan)

    ! Normalize the window function.
    normalization=windowFunction1(1,1,1)
    if (real(normalization) > 0.0d0) windowFunction1=windowFunction1/normalization
    normalization=windowFunction2(1,1,1)
    if (real(normalization) > 0.0d0) windowFunction2=windowFunction2/normalization

    return
  end subroutine Geometry_Survey_Window_Functions_Li_White_2009_SDSS
  
end module Geometry_Surveys_Li_White_2009_SDSS
