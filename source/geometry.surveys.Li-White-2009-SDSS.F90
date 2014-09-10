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

  type, extends(surveyGeometryRandomPoints) :: surveyGeometryLiWhite2009SDSS
   contains
     procedure :: distanceMaximum   => liWhite2009SDSSDistanceMaximum
     procedure :: solidAngle        => liWhite2009SDSSSolidAngle
     procedure :: randomsInitialize => liWhite2009SDSSRandomsInitialize
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

    liWhite2009SDSSDefaultConstructor%geometryInitialized=.false.
    return
  end function liWhite2009SDSSDefaultConstructor

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

  subroutine liWhite2009SDSSRandomsInitialize(self)
    !% Compute the window function for the survey.
    use ISO_Varying_String
    use String_Handling
    use File_Utilities
    use Galacticus_Input_Paths
    use Memory_Management
    use System_Command
    use Galacticus_Error
    use Galacticus_Display
    use Numerical_Constants_Math
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)             :: self
    double precision                               , allocatable, dimension(:) :: angleTmp
    integer                                                                    :: randomsCount  , j          , &
         &                                                                        i             , randomUnit
    double precision                                                           :: rightAscension, declination
    type            (varying_string               )                            :: message
    
    ! Randoms file obtained from:  http://sdss.physics.nyu.edu/lss/dr72/random/
    if (.not.File_Exists(Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")) then
       call System_Command_Do("mkdir -p "//Galacticus_Input_Path()//"data/surveyGeometry")
       call System_Command_Do("wget http://sdss.physics.nyu.edu/lss/dr72/random/lss_random-0.dr72.dat -O "//Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")
       if (.not.File_Exists(Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")) call Galacticus_Error_Report('liWhite2009SDSSWindowFunctions','unable to download SDSS survey geometry randoms file')
    end if
    randomsCount=Count_Lines_In_File(Galacticus_Input_Path()//"data/surveyGeometry/lss_random-0.dr72.dat")
    call Alloc_Array(self%randomTheta,[randomsCount])
    call Alloc_Array(self%randomPhi  ,[randomsCount])
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
          self%randomTheta(j)=Pi*(90.0d0-declination   )/180.0d0
          self%randomPhi  (j)=Pi*        rightAscension /180.0d0
       end if
    end do
    close(randomUnit)
    randomsCount=j
    call Move_Alloc   (self%randomTheta,angleTmp      )
    call Alloc_Array  (self%randomTheta,[randomsCount])
    self%randomTheta=angleTmp(1:randomsCount)
    call Dealloc_Array(angleTmp                       )
    call Move_Alloc   (self%randomPhi  ,angleTmp      )
    call Alloc_Array  (self%randomPhi  ,[randomsCount])
    self%randomPhi  =angleTmp(1:randomsCount)
    call Dealloc_Array(angleTmp                       )
    message="Read "
    message=message//randomsCount//" random points and kept "//randomsCount//" of them"
    call Galacticus_Display_Message(message)
    return
  end subroutine liWhite2009SDSSRandomsInitialize
