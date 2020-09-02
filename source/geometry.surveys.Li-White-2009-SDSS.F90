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

!% Implements the survey geometry of the SDSS sample used by \cite{li_distribution_2009}.

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !# <surveyGeometry name="surveyGeometryLiWhite2009SDSS">
  !#  <description>Implements the survey geometry of the SDSS sample used by \cite{li_distribution_2009}.</description>
  !# </surveyGeometry>
  type, extends(surveyGeometryRandomPoints) :: surveyGeometryLiWhite2009SDSS
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_  => null()
     double precision                                   :: redshiftMinimum               , redshiftMaximum
     double precision                                   :: limitDistanceMinimum          , limitDistanceMaximum
   contains
     final     ::                      liWhite2009SDSSDestructor
     procedure :: distanceMinimum   => liWhite2009SDSSDistanceMinimum
     procedure :: distanceMaximum   => liWhite2009SDSSDistanceMaximum
     procedure :: solidAngle        => liWhite2009SDSSSolidAngle
     procedure :: randomsInitialize => liWhite2009SDSSRandomsInitialize
  end type surveyGeometryLiWhite2009SDSS

  interface surveyGeometryLiWhite2009SDSS
     !% Constructors for the \cite{li_distribution_2009} survey geometry class.
     module procedure liWhite2009SDSSConstructorParameters
     module procedure liWhite2009SDSSConstructorInternal
  end interface surveyGeometryLiWhite2009SDSS

contains

  function liWhite2009SDSSConstructorParameters(parameters) result(self)
    !% Constructor for the \cite{li_distribution_2009} survey geometry class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (surveyGeometryLiWhite2009SDSS)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (randomNumberGeneratorClass   ), pointer       :: randomNumberGenerator_
    double precision                                               :: redshiftMinimum       , redshiftMaximum

    !# <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !# <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !# <inputParameter>
    !#   <name>redshiftMinimum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The minimum redshift for the survey.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshiftMaximum</name>
    !#   <defaultValue>huge(1.0d0)</defaultValue>
    !#   <source>parameters</source>
    !#   <description>The maximum redshift for the survey.</description>
    !# </inputParameter>
    ! Build the object.
    self=surveyGeometryLiWhite2009SDSS(redshiftMinimum,redshiftMaximum,cosmologyFunctions_,randomNumberGenerator_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"   />
    !# <objectDestructor name="randomNumberGenerator_"/>
    return
  end function liWhite2009SDSSConstructorParameters

  function liWhite2009SDSSConstructorInternal(redshiftMinimum,redshiftMaximum,cosmologyFunctions_,randomNumberGenerator_) result(self)
    !% Constructor for the \cite{li_distribution_2009} survey geometry class which allows specification of minimum and maximum redshifts.
    use :: Cosmology_Functions        , only : cosmologyFunctionsClass
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    implicit none
    type            (surveyGeometryLiWhite2009SDSS)                                  :: self
    double precision                               , intent(in   )                   :: redshiftMinimum    , redshiftMaximum
    class           (cosmologyFunctionsClass      ), intent(in   ), target           :: cosmologyFunctions_
    class           (randomNumberGeneratorClass   ), intent(in   ), target, optional :: randomNumberGenerator_ 
    !# <constructorAssign variables="*cosmologyFunctions_, *randomNumberGenerator_, redshiftMinimum, redshiftMaximum"/>

    self   %geometryInitialized =.false.
    self   %limitDistanceMinimum=self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                                        output  =distanceTypeComoving, &
         &                                                                        redshift=redshiftMinimum       &
         &                                                                       )
    if (redshiftMaximum < huge(1.0d0)) then
       self%limitDistanceMaximum=self%cosmologyFunctions_%distanceComovingConvert(                               &
            &                                                                     output  =distanceTypeComoving, &
            &                                                                     redshift=redshiftMaximum       &
            &                                                                    )
    else
       self%limitDistanceMaximum=huge(1.0d0)
    end if
    return
  end function liWhite2009SDSSConstructorInternal

  subroutine liWhite2009SDSSDestructor(self)
    !% Destructor for the ``liWhite2009SDSS'' survey geometry class.
    implicit none
    type(surveyGeometryLiWhite2009SDSS), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%randomNumberGenerator_"/>
    return
  end subroutine liWhite2009SDSSDestructor

  double precision function liWhite2009SDSSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,field)
    !% Compute the minimum distance at which a galaxy is visible.
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)           :: self
    double precision                               , intent(in   ), optional :: mass , magnitudeAbsolute, luminosity
    integer                                        , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity

    liWhite2009SDSSDistanceMinimum=self%limitDistanceMinimum
    return
  end function liWhite2009SDSSDistanceMinimum

  double precision function liWhite2009SDSSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Galacticus_Error           , only : Galacticus_Error_Report
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)           :: self
    double precision                               , intent(in   ), optional :: mass    , magnitudeAbsolute, luminosity
    integer                                        , intent(in   ), optional :: field
    double precision                                                         :: redshift, logarithmicMass
    !$GLC attributes unused :: self, magnitudeAbsolute, luminosity

    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('field = 1 required'//{introspection:location})
    ! Find the limiting redshift for this mass using a fit derived from Millennium Simulation SAMs. (See
    ! constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07/massLuminosityRelation.pl for details.)
    logarithmicMass=log10(mass)
    redshift=                                   &
         & max(                                 &
         &     -5.9502006195004d0               &
         &     +logarithmicMass                 &
         &     *(                               &
         &       +2.63793788603951d0            &
         &       +logarithmicMass               &
         &       *(                             &
         &         -0.421075858899237d0         &
         &         +logarithmicMass             &
         &         *(                           &
         &           +0.0285198776926787d0      &
         &           +logarithmicMass           &
         &           *(                         &
         &             -0.000678327494720407d0  &
         &            )                         &
         &          )                           &
         &        )                             &
         &      )                             , &
         &     +0.0d0                           &
         &    )
    ! Convert from redshift to comoving distance.
    liWhite2009SDSSDistanceMaximum=min(                                                                                &
         &                             self%limitDistanceMaximum                                                     , &
         &                             self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                                              output  =distanceTypeComoving, &
         &                                                                              redshift=redshift              &
         &                                                                             )                               &
         &                            )
    return
  end function liWhite2009SDSSDistanceMaximum

  double precision function liWhite2009SDSSSolidAngle(self,field)
    !% Return the solid angle of the \cite{li_distribution_2009} sample.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)           :: self
    integer                                        , intent(in   ), optional :: field
    double precision                               , parameter               :: solidAngleSurvey=2.1901993d0 ! From Percival et al. (2010; MNRAS; 401; 2148)
    !$GLC attributes unused :: self

    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('field = 1 required'//{introspection:location})
    liWhite2009SDSSSolidAngle=solidAngleSurvey
    return
  end function liWhite2009SDSSSolidAngle

  subroutine liWhite2009SDSSRandomsInitialize(self)
    !% Compute the window function for the survey.
    use :: File_Utilities          , only : Count_Lines_In_File       , Directory_Make     , File_Exists
    use :: Galacticus_Display      , only : Galacticus_Display_Message
    use :: Galacticus_Error        , only : Galacticus_Error_Report
    use :: Galacticus_Paths        , only : galacticusPath            , pathTypeDataDynamic
    use :: ISO_Varying_String      , only : varying_string
    use :: Memory_Management       , only : allocateArray             , deallocateArray
    use :: Numerical_Constants_Math, only : Pi
    use :: String_Handling         , only : operator(//)
    use :: System_Command          , only : System_Command_Do
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)             :: self
    double precision                               , allocatable, dimension(:) :: angleTmp
    integer                                                                    :: randomsCount  , j          , &
         &                                                                        i             , randomUnit
    double precision                                                           :: rightAscension, declination
    type            (varying_string               )                            :: message

    ! Randoms file obtained from:  http://sdss.physics.nyu.edu/lss/dr72/random/
    if (.not.File_Exists(galacticusPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat")) then
       call Directory_Make(galacticusPath(pathTypeDataDynamic)//"surveyGeometry")
       call System_Command_Do("wget http://sdss.physics.nyu.edu/lss/dr72/random/lss_random-0.dr72.dat -O "//galacticusPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat")
       if (.not.File_Exists(galacticusPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat")) call Galacticus_Error_Report('unable to download SDSS survey geometry randoms file'//{introspection:location})
    end if
    randomsCount=Count_Lines_In_File(galacticusPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat")
    call allocateArray(self%randomTheta,[randomsCount])
    call allocateArray(self%randomPhi  ,[randomsCount])
    open(newUnit=randomUnit,file=char(galacticusPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat"),status="old",form="formatted")
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
    call allocateArray  (self%randomTheta,[randomsCount])
    self%randomTheta=angleTmp(1:randomsCount)
    call deallocateArray(angleTmp                       )
    call Move_Alloc   (self%randomPhi  ,angleTmp      )
    call allocateArray  (self%randomPhi  ,[randomsCount])
    self%randomPhi  =angleTmp(1:randomsCount)
    call deallocateArray(angleTmp                       )
    message="Read "
    message=message//randomsCount//" random points and kept "//randomsCount//" of them"
    call Galacticus_Display_Message(message)
    return
  end subroutine liWhite2009SDSSRandomsInitialize
