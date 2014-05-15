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

!% Implements the survey geometry of the SDSS sample used by \cite{li_distribution_2009}.

  !# <surveyGeometry name="surveyGeometryLiWhite2009SDSS">
  !#  <description>Implements the survey geometry of the SDSS sample used by \cite{li_distribution_2009}.</description>
  !# </surveyGeometry>

  type, extends(surveyGeometryClass) :: surveyGeometryLiWhite2009SDSS
   contains
     procedure :: windowFunctionAvailable => liWhite2009SDSSWindowFunctionAvailable
     procedure :: angularPowerAvailable   => liWhite2009SDSSAngularPowerAvailable
     procedure :: distanceMaximum         => liWhite2009SDSSDistanceMaximum
     procedure :: solidAngle              => liWhite2009SDSSSolidAngle
     procedure :: windowFunctions         => liWhite2009SDSSWindowFunctions
     procedure :: angularPower            => liWhite2009SDSSAngularPower
  end type surveyGeometryLiWhite2009SDSS

  interface surveyGeometryLiWhite2009SDSS
     !% Constructors for the \cite{li_distribution_2009} survey geometry class.
     module procedure liWhite2009SDSSDefaultConstructor
  end interface surveyGeometryLiWhite2009SDSS

contains

  function liWhite2009SDSSDefaultConstructor()
    !% Default constructor for the \cite{li_distribution_2009} conditional mass function class.
    use Input_Parameters
    implicit none
    type(surveyGeometryLiWhite2009SDSS) :: liWhite2009SDSSDefaultConstructor

    return
  end function liWhite2009SDSSDefaultConstructor

  logical function liWhite2009SDSSWindowFunctionAvailable(self)
    !% Return true to indicate that survey window function is available.
    implicit none
    class(surveyGeometryLiWhite2009SDSS), intent(inout) :: self

    liWhite2009SDSSWindowFunctionAvailable=.true.
    return
  end function liWhite2009SDSSWindowFunctionAvailable

  logical function liWhite2009SDSSAngularPowerAvailable(self)
    !% Return false to indicate that survey angular power is not available.
    implicit none
    class(surveyGeometryLiWhite2009SDSS), intent(inout) :: self

    liWhite2009SDSSAngularPowerAvailable=.false.
    return
  end function liWhite2009SDSSAngularPowerAvailable

  double precision function liWhite2009SDSSDistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    use Galacticus_Error
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)           :: self
    double precision                               , intent(in   )           :: mass
    integer                                        , intent(in   ), optional :: field
    class           (cosmologyFunctionsClass      ), pointer                 :: cosmologyFunctions_
    double precision                                                         :: redshift           , logarithmicMass
    
!! AJB HACK
type(surveyGeometryBernardi2013SDSS) :: b
b=surveyGeometryBernardi2013SDSS()
liWhite2009SDSSDistanceMaximum=b%distanceMaximum(mass,field)
return


    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('liWhite2009SDSSDistanceMaximum','field = 1 required')
    ! Find the limiting redshift for this mass using a fit derived from Millennium Simulation SAMs. (See
    ! constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07/massLuminosityRelation.pl for details.)
    logarithmicMass=log10(mass)
    redshift=                              &
         & -5.9502006195004d0              &
         & +logarithmicMass                &
         & *(                              &
         &   +2.63793788603951d0           &
         &   +logarithmicMass              &
         &   *(                            &
         &     -0.421075858899237d0        &
         &     +logarithmicMass            &
         &     *(                          &
         &       +0.0285198776926787d0     &
         &       +logarithmicMass          &
         &       *(                        &
         &         -0.000678327494720407d0 &
         &        )                        &
         &      )                          &
         &    )                            &
         &  )
    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()     
    ! Convert from redshift to comoving distance.
    liWhite2009SDSSDistanceMaximum                                                     &
         & =cosmologyFunctions_%distanceComovingConvert(                               &
         &                                              output  =distanceTypeComoving, &
         &                                              redshift=redshift              &
         &                                             )
    return
  end function liWhite2009SDSSDistanceMaximum

  double precision function liWhite2009SDSSSolidAngle(self,field)
    !% Return the solid angle of the \cite{li_distribution_2009} sample.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)           :: self
    integer                                        , intent(in   ), optional :: field
    double precision                               , parameter               :: solidAngleSurvey=2.1901993d0 ! From Percival et al. (2010; MNRAS; 401; 2148)
    
    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('liWhite2009SDSSSolidAngle','field = 1 required')
    liWhite2009SDSSSolidAngle=solidAngleSurvey
    return
  end function liWhite2009SDSSSolidAngle

  subroutine liWhite2009SDSSWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
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
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)                                               :: self
    double precision                               , intent(in   )                                               :: mass1,mass2
    integer                                        , intent(in   )                                               :: gridCount
    double precision                               , intent(  out)                                               :: boxLength
    complex         (c_double_complex             ), intent(  out),     dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2
    double precision                               ,                    dimension(3                            ) :: origin,position1,position2
    integer                                                                                                      :: i,j
    double precision                                                                                             :: comovingDistanceMaximum1, comovingDistanceMaximum2, &
         &                                                                                                          comovingDistanceMinimum1, comovingDistanceMinimum2
    type            (c_ptr                        )                                                              :: plan
    complex         (c_double_complex             ),                    dimension(gridCount,gridCount,gridCount) :: selectionFunction1,selectionFunction2
    complex         (c_double_complex             )                                                              :: normalization
    logical                                        , save                                                        :: geometryInitialized=.false.
    double precision                                                                                             :: rightAscension,declination,distance1,distance2
    double precision                               , save, allocatable, dimension(:                            ) :: randomTheta,randomPhi
    integer                                        , save                                                        :: randomsCount
    integer                                                                                                      :: randomUnit
    type            (fgsl_rng                     ), save                                                        :: pseudoSequenceObject
    logical                                        , save                                                        :: reset=.true.
    type            (varying_string               )                                                              :: message
    double precision                               , save                                                        :: surveyDistanceMinimum&
         &,surveyDistanceMaximum
    class           (cosmologyFunctionsClass      ), pointer                                                     :: cosmologyFunctions_

    ! Initialize geometry if necessary.
    if (.not.geometryInitialized) then
       ! Randoms file obtained from:  http://sdss.physics.nyu.edu/lss/dr72/random/
       if (.not.File_Exists(Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")) then
          call System_Command_Do("mkdir -p "//Galacticus_Input_Path()//"data/surveyGeometry")
          call System_Command_Do("wget http://sdss.physics.nyu.edu/lss/dr72/random/lss_random-0.dr72.dat -O "//Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")
          if (.not.File_Exists(Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")) call Galacticus_Error_Report('liWhite2009SDSSWindowFunctions','unable to download SDSS survey geometry randoms file')
       end if
       randomsCount=Count_Lines_In_File(Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")
       call Alloc_Array(randomTheta,[randomsCount])
       call Alloc_Array(randomPhi  ,[randomsCount])
       open(newUnit=randomUnit,file=char(Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat"),status="old",form="formatted")
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

       ! Get the default cosmology functions object.
       cosmologyFunctions_ => cosmologyFunctions()    
       ! Compute the distances corresponding to the minimum and maximum redshifts.
       surveyDistanceMinimum=cosmologyFunctions_%distanceComoving(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(0.001d0)))
       surveyDistanceMaximum=cosmologyFunctions_%distanceComoving(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(0.500d0)))
       geometryInitialized=.true.
    end if

    ! Find the comoving distance corresponding to this distance module.
    comovingDistanceMaximum1=min(self%distanceMaximum(mass1),surveyDistanceMaximum)
    comovingDistanceMaximum2=min(self%distanceMaximum(mass2),surveyDistanceMaximum)
    comovingDistanceMinimum1=surveyDistanceMinimum
    comovingDistanceMinimum2=surveyDistanceMinimum

    ! Determine a suitable box length for the window function calculation.
    boxLength=3.0d0*max(comovingDistanceMaximum1,comovingDistanceMaximum2)

    ! Set up origin.
    origin=([0.0d0,0.5d0,0.5d0]/dble(gridCount)+[0.5d0,0.5d0,0.5d0])*dble(boxLength)

    ! Populate the cube with the survey selection function.
    selectionFunction1=cmplx(0.0d0,0.0d0,kind=c_double_complex)
    selectionFunction2=cmplx(0.0d0,0.0d0,kind=c_double_complex)

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
  end subroutine liWhite2009SDSSWindowFunctions

  double precision function liWhite2009SDSSAngularPower(self,i,j,l)
    !% Angular power is not available, so simply aborts.
    use Galacticus_Error
    implicit none
    class  (surveyGeometryLiWhite2009SDSS), intent(inout) :: self
    integer                               , intent(in   ) :: i   , j, l

    call Galacticus_Error_Report('liWhite2009SDSSAngularPower','angular power is not available')
    return
  end function liWhite2009SDSSAngularPower
