!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+    Contributions to this file made by:  Luiz Felippe S. Rodrigues.

  !% An implementation of the intergalactic medium state class in which state is read from file.

  !# <intergalacticMediumState name="intergalacticMediumStateFile">
  !#  <description>The intergalactic medium state is read from file.</description>
  !# </intergalacticMediumState>
  use FGSL

  ! Current file format version for intergalactic medium state files.
  integer, parameter :: fileFormatVersionCurrent=1

  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateFile
     !% An \gls{igm} state class which reads state from file.
     private
     ! Name of the file from which to read intergalactic medium state data.
     type            (varying_string   )                            :: fileName     
     ! Flag indicating whether or not data has been read.
     logical                                                        :: dataRead                                       =.false.
     ! Data tables.
     integer                                                        :: redshiftCount
     double precision                   , allocatable, dimension(:) :: electronFractionTable                                  , temperatureTable                                    , &
          &                                                            timeTable                                              , ionizedHydrogenFractionTable                        , &
          &                                                            ionizedHeliumFractionTable
     ! Interpolation objects.
     type            (fgsl_interp_accel)                            :: interpolationAcceleratorElectronFraction               , interpolationAcceleratorTemperature                 , &
          &                                                            interpolationAcceleratorIonizedHydrogenFraction        , interpolationAcceleratorIonizedHeliumFraction
     type            (fgsl_interp      )                            :: interpolationObjectElectronFraction                    , interpolationObjectTemperature                      , &
          &                                                            interpolationObjectIonizedHydrogenFraction             , interpolationObjectIonizedHeliumFraction
     logical                                                        :: interpolationResetElectronFraction             =.true. , interpolationResetTemperature                =.true., &
          &                                                            interpolationResetIonizedHydrogenFraction      =.true. , interpolationResetIonizedHeliumFraction      =.true.
   contains
     procedure :: electronFraction            => fileElectronFraction
     procedure :: neutralHydrogenFraction     => fileNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => fileNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => fileSinglyIonizedHeliumFraction
     procedure :: temperature                 => fileTemperature
  end type intergalacticMediumStateFile
  
  interface intergalacticMediumStateFile
     !% Constructors for the file intergalactic medium state class.
     module procedure fileDefaultConstructor
     module procedure fileConstructor
  end interface intergalacticMediumStateFile

contains

  function fileDefaultConstructor()
    !% Default constructor for the file \gls{igm} state class.
    use Input_Parameters
    implicit none
    type(intergalacticMediumStateFile), target  :: fileDefaultConstructor
    type(varying_string              )          :: intergalaticMediumStateFileName
    
    !@ <inputParameter>
    !@   <name>intergalaticMediumStateFileName</name>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The name of the file from which to read intergalactic medium state data.
    !@   </description>
    !@   <type>string</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('intergalaticMediumStateFileName',intergalaticMediumStateFileName)
    ! Construct the object.
    fileDefaultConstructor=fileConstructor(intergalaticMediumStateFileName)
    return
  end function fileDefaultConstructor

  function fileConstructor(fileName)
    !% Constructor for the file \gls{igm} state class.
    use Cosmology_Functions
    implicit none
    type(intergalacticMediumStateFile), target        :: fileConstructor
    type(varying_string              ), intent(in   ) :: fileName

    fileConstructor%fileName=fileName
    return
  end function fileConstructor
  
  subroutine fileReadData(self)
    !% Read in data describing the state of the intergalactic medium.
    use Galacticus_Error
    use FoX_dom
    use IO_XML
    use Cosmology_Functions
    implicit none
    class  (intergalacticMediumStateFile), intent(inout) :: self
    type   (node                        ), pointer       :: doc                      , thisItem
    class  (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctionsDefault
    integer                                              :: fileFormatVersion        , iRedshift, &
         &                                                  ioStatus

    ! Check if data has yet to be read.
    if (.not.self%dataRead) then
       !$omp critical (FoX_DOM_Access)
       doc => parseFile(char(self%fileName),iostat=ioStatus)
       if (ioStatus /= 0) call Galacticus_Error_Report('fileReadData','Unable to parse intergalactic medium state file')
       ! Check the file format version of the file.
       thisItem             => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
       call extractDataContent(thisItem,fileFormatVersion)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('fileReadData','file format version is out of date')
       ! Read the data.
       thisItem             => XML_Get_First_Element_By_Tag_Name(doc,"redshift"          )
       call XML_Array_Read(thisItem,"datum",self%timeTable                   )
       thisItem             => XML_Get_First_Element_By_Tag_Name(doc,"electronFraction"  )
       call XML_Array_Read(thisItem,"datum",self%electronFractionTable       )
       thisItem             => XML_Get_First_Element_By_Tag_Name(doc,"hIonizedFraction"  )
       call XML_Array_Read(thisItem,"datum",self%ionizedHydrogenFractionTable)
       thisItem             => XML_Get_First_Element_By_Tag_Name(doc,"heIonizedFraction" )
       call XML_Array_Read(thisItem,"datum",self%ionizedHeliumFractionTable  )
       thisItem             => XML_Get_First_Element_By_Tag_Name(doc,"matterTemperature" )
       call XML_Array_Read(thisItem,"datum",self%temperatureTable     )
       self%redshiftCount=size(self%timeTable)
       ! Get the default cosmology functions object.
       cosmologyFunctionsDefault => cosmologyFunctions()
       ! Convert redshifts to times.
       do iRedshift=1,self%redshiftCount
          self%timeTable(iRedshift)=cosmologyFunctionsDefault%cosmicTime(cosmologyFunctionsDefault%expansionFactorFromRedshift(self%timeTable(iRedshift)))
       end do
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
       ! Flag that data has now been read.
       self%dataRead=.true.
    end if
    return
  end subroutine fileReadData

  double precision function fileTemperature(self,time)
    !% Return the temperature of the intergalactic medium at the specified {\tt time} by interpolating in tabulated data.
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the electron fraction.
    fileTemperature                                               &
         & =Interpolate(                                          &
         &              self%redshiftCount                      , &
         &              self%timeTable                          , &
         &              self%temperatureTable                   , &
         &              self%interpolationObjectTemperature     , &
         &              self%interpolationAcceleratorTemperature, &
         &              time                                    , &
         &              reset=self%interpolationResetTemperature  &
         &             )
    return
  end function fileTemperature

  double precision function fileElectronFraction(self,time)
    !% Return the electron fraction in the intergalactic medium at the specified {\tt time} by interpolating in tabulated data,
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the electron fraction.
    fileElectronFraction                                               &
         & =Interpolate(                                               &
         &              self%redshiftCount                           , &
         &              self%timeTable                               , &
         &              self%electronFractionTable                   , &
         &              self%interpolationObjectElectronFraction     , &
         &              self%interpolationAcceleratorElectronFraction, &
         &              time                                         , &
         &              reset=self%interpolationResetElectronFraction  &
         &             )
    return
  end function fileElectronFraction

  double precision function fileNeutralHydrogenFraction(self,time)
    !% Return the neutral hydrogen fraction in the intergalactic medium at the specified {\tt time} by interpolating in tabulated data,
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the electron fraction.
    fileNeutralHydrogenFraction                                                &
         & =+1.0d0                                                             &
         &  -Interpolate(                                                      &
         &               self%redshiftCount                                  , &
         &               self%timeTable                                      , &
         &               self%ionizedHydrogenFractionTable                   , &
         &               self%interpolationObjectIonizedHydrogenFraction     , &
         &               self%interpolationAcceleratorIonizedHydrogenFraction, &
         &               time                                                , &
         &               reset=self%interpolationResetIonizedHydrogenFraction  &
         &              )
    return
  end function fileNeutralHydrogenFraction

  double precision function fileNeutralHeliumFraction(self,time)
    !% Return the neutral helium fraction in the intergalactic medium at the specified {\tt time} by interpolating in tabulated data,
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the electron fraction.
    fileNeutralHeliumFraction                                                &
         & =+1.0d0                                                           &
         &  -Interpolate(                                                    &
         &               self%redshiftCount                                , &
         &               self%timeTable                                    , &
         &               self%ionizedHeliumFractionTable                   , &
         &               self%interpolationObjectIonizedHeliumFraction     , &
         &               self%interpolationAcceleratorIonizedHeliumFraction, &
         &               time                                              , &
         &               reset=self%interpolationResetIonizedHeliumFraction  &
         &              )
    return
  end function fileNeutralHeliumFraction

  double precision function fileSinglyIonizedHeliumFraction(self,time)
    !% Return the neutral helium fraction in the intergalactic medium at the specified {\tt time} by interpolating in tabulated data,
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the electron fraction.
    fileSinglyIonizedHeliumFraction                                         &
         & =Interpolate(                                                    &
         &              self%redshiftCount                                , &
         &              self%timeTable                                    , &
         &              self%ionizedHeliumFractionTable                   , &
         &              self%interpolationObjectIonizedHeliumFraction     , &
         &              self%interpolationAcceleratorIonizedHeliumFraction, &
         &              time                                              , &
         &              reset=self%interpolationResetIonizedHeliumFraction  &
         &             )
    return
  end function fileSinglyIonizedHeliumFraction
