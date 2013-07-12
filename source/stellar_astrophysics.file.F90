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

!% Contains a module which implements calculation related to stellar astrophyics.

module Stellar_Astrophysics_File
  !% Implements calculation related to stellar astrophyics.
  use Numerical_Interpolation_2D_Irregular
  implicit none
  private
  public :: Stellar_Astrophysics_File_Initialize, Stellar_Astrophysics_File_Format_Version

  ! Arrays to store stellar properties.
  double precision, allocatable, dimension(:  ) :: stellarLifetime             , stellarLifetimeMass, &
       &                                           stellarLifetimeMetallicity
  double precision, allocatable, dimension(:  ) :: ejectedMass                 , ejectedMassMass    , &
       &                                           ejectedMassMetallicity
  double precision, allocatable, dimension(:  ) :: metalYield                  , metalYieldMass     , &
       &                                           metalYieldMetallicity
  double precision, allocatable, dimension(:,:) :: elementYield                , elementYieldMass   , &
       &                                           elementYieldMetallicity

  ! Variables that store information about number of elements and number of yields for each element.
  integer         , allocatable, dimension(:)   :: atomIndexMap                , elementYieldCount

  ! Number of elements being tracked.
  integer                                       :: elementCount

  ! Current file format version for intergalactic background radiation files.
  integer         , parameter                   :: fileFormatVersionCurrent  =1

contains

  integer function Stellar_Astrophysics_File_Format_Version()
    !% Return the current file format version of stellar astrophysics files.
    implicit none

    Stellar_Astrophysics_File_Format_Version=fileFormatVersionCurrent
    return
  end function Stellar_Astrophysics_File_Format_Version

  !# <stellarAstrophysicsMethod>
  !#  <unitName>Stellar_Astrophysics_File_Initialize</unitName>
  !# </stellarAstrophysicsMethod>
  subroutine Stellar_Astrophysics_File_Initialize(stellarAstrophysicsMethod,Star_Ejected_Mass_Get,Star_Initial_Mass_Get&
       &,Star_Metal_Yield_Mass_Get,Star_Lifetime_Get)
    !% Initialize the stellar astrophysics module.
    use FoX_dom
    use Galacticus_Error
    use Memory_Management
    use ISO_Varying_String
    use Input_Parameters
    use Atomic_Data
    use Galacticus_Input_Paths
    implicit none
    type            (varying_string            ), intent(in   )          :: stellarAstrophysicsMethod
    procedure       (Star_Ejected_Mass_File    ), intent(inout), pointer :: Star_Ejected_Mass_Get
    procedure       (Star_Initial_Mass_File    ), intent(inout), pointer :: Star_Initial_Mass_Get
    procedure       (Star_Metal_Yield_Mass_File), intent(inout), pointer :: Star_Metal_Yield_Mass_Get
    procedure       (Star_Lifetime_File        ), intent(inout), pointer :: Star_Lifetime_Get
    type            (Node                      )               , pointer :: doc                      , thisDatum               , &
         &                                                                  thisStar
    type            (NodeList                  )               , pointer :: propertyList             , starList
    type            (varying_string            )                         :: stellarPropertiesFile
    integer                                                              :: ejectedMassCount         , elementYieldCountMaximum, &
         &                                                                  fileFormatVersion        , iElement                , &
         &                                                                  iStar                    , ioErr                   , &
         &                                                                  lifetimeCount            , mapToIndex              , &
         &                                                                  metalYieldCount
    double precision                                                     :: initialMass              , metallicity
    logical                                                              :: starHasElements

    ! Check if our method is selected.
    if (stellarAstrophysicsMethod == 'file') then
       ! Set up procedure pointers.
       Star_Ejected_Mass_Get     => Star_Ejected_Mass_File
       Star_Initial_Mass_Get     => Star_Initial_Mass_File
       Star_Metal_Yield_Mass_Get => Star_Metal_Yield_Mass_File
       Star_Lifetime_Get         => Star_Lifetime_File

       ! Get the name of the file containing stellar data.
       !@ <inputParameter>
       !@   <name>stellarPropertiesFile</name>
       !@   <defaultValue>data/stellarAstrophysics/Stellar\_Properties\_Compilation.xml</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the XML file from which to read stellar properties (ejected masses, yields, etc.).
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPropertiesFile',stellarPropertiesFile,defaultValue=char(Galacticus_Input_Path())//'data&
            &/stellarAstrophysics/Stellar_Properties_Compilation.xml')

       ! Allocate array to store number of entries in file for yield of each element.
       call Alloc_Array(elementYieldCount,[Atomic_Data_Atoms_Count()])
       call Alloc_Array(atomIndexMap     ,[Atomic_Data_Atoms_Count()])

       !$omp critical (FoX_DOM_Access)
       ! Open the XML file containing stellar properties.
       doc => parseFile(char(stellarPropertiesFile),iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize','Unable to parse stellar properties file')

       ! Check the file format version of the file.
       propertyList => getElementsByTagname(doc,"fileFormat")
       thisDatum => item(propertyList,0)
       call extractDataContent(thisDatum,fileFormatVersion)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('Stellar_Astrophysics_File_Initialize','file format version is out of date')

       ! Get a list of all stars.
       starList => getElementsByTagname(doc,"star")

       ! Count up number of stars with given properties.
       lifetimeCount    =0
       ejectedMassCount =0
       metalYieldCount  =0
       elementYieldCount=0
       do iStar=0,getLength(starList)-1
          thisStar     => item(starList,iStar)
          propertyList => getElementsByTagname(thisStar,"initialMass")
          if (getLength(propertyList) /= 1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize'&
               & ,'star must have precisely one initial mass')
          propertyList => getElementsByTagname(thisStar,"metallicity")
          if (getLength(propertyList) /= 1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize'&
               & ,'star must have precisely one metallicity')
          propertyList => getElementsByTagname(thisStar,"lifetime")
          if (getLength(propertyList) == 1) lifetimeCount=lifetimeCount+1
          if (getLength(propertyList) >  1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize'&
               & ,'star has multiple lifetimes')
          propertyList => getElementsByTagname(thisStar,"ejectedMass")
          if (getLength(propertyList) == 1) ejectedMassCount=ejectedMassCount+1
          if (getLength(propertyList) >  1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize'&
               & ,'star has multiple ejected masses')
          propertyList => getElementsByTagname(thisStar,"metalYieldMass")
          if (getLength(propertyList) == 1) metalYieldCount=metalYieldCount+1
          if (getLength(propertyList) >  1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize'&
               & ,'star has multiple metal yield masses')
          do iElement=1,size(elementYieldCount)
             propertyList => getElementsByTagname(thisStar,"elementYieldMass"//trim(Atomic_Short_Label(iElement)))
             if (getLength(propertyList) == 1) elementYieldCount(iElement)=elementYieldCount(iElement)+1
             if (getLength(propertyList) >  1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize'&
                  & ,'star has multiple element yield masses')
          end do
       end do

       ! Find number of elements for which some yield data is available.
       elementCount            =count (elementYieldCount > 0)
       elementYieldCountMaximum=maxval(elementYieldCount    )

       ! Create mapping of atomic index to our array space.
       mapToIndex=0
       do iElement=1,size(elementYieldCount)
          if (elementYieldCount(iElement) > 0) then
             mapToIndex=mapToIndex+1
             atomIndexMap(iElement)=mapToIndex
          else
             atomIndexMap(iElement)=-1
          end if
       end do

       ! Allocate arrays to store stellar properties.
       call Alloc_Array(stellarLifetime           ,[lifetimeCount                        ])
       call Alloc_Array(stellarLifetimeMass       ,[lifetimeCount                        ])
       call Alloc_Array(stellarLifetimeMetallicity,[lifetimeCount                        ])
       call Alloc_Array(ejectedMass               ,[ejectedMassCount                     ])
       call Alloc_Array(ejectedMassMass           ,[ejectedMassCount                     ])
       call Alloc_Array(ejectedMassMetallicity    ,[ejectedMassCount                     ])
       call Alloc_Array(metalYield                ,[metalYieldCount                      ])
       call Alloc_Array(metalYieldMass            ,[metalYieldCount                      ])
       call Alloc_Array(metalYieldMetallicity     ,[metalYieldCount                      ])
       call Alloc_Array(elementYield              ,[elementYieldCountMaximum,elementCount])
       call Alloc_Array(elementYieldMass          ,[elementYieldCountMaximum,elementCount])
       call Alloc_Array(elementYieldMetallicity   ,[elementYieldCountMaximum,elementCount])

       ! Loop over stars to process their properties.
       lifetimeCount    =0
       ejectedMassCount =0
       metalYieldCount  =0
       elementYieldCount=0
       do iStar=0,getLength(starList)-1
          thisStar     => item(starList,iStar)
          propertyList => getElementsByTagname(thisStar,"initialMass")
          thisDatum    => item(propertyList,0)
          call extractDataContent(thisDatum,initialMass)
          propertyList => getElementsByTagname(thisStar,"metallicity")
          thisDatum    => item(propertyList,0)
          call extractDataContent(thisDatum,metallicity)

          ! Process stellar lifetimes.
          propertyList => getElementsByTagname(thisStar,"lifetime")
          if (getLength(propertyList) == 1) then
             thisDatum => item(propertyList,0)
             lifetimeCount=lifetimeCount+1
             call extractDataContent(thisDatum,stellarLifetime(lifetimeCount))
             stellarLifetimeMass       (lifetimeCount)=initialMass
             stellarLifetimeMetallicity(lifetimeCount)=metallicity
          end if

          ! Process ejected masses.
          propertyList => getElementsByTagname(thisStar,"ejectedMass")
          if (getLength(propertyList) == 1) then
             thisDatum => item(propertyList,0)
             ejectedMassCount=ejectedMassCount+1
             call extractDataContent(thisDatum,ejectedMass(ejectedMassCount))
             ejectedMassMass       (ejectedMassCount)=initialMass
             ejectedMassMetallicity(ejectedMassCount)=metallicity
          end if

          ! Process metal yields.
          propertyList => getElementsByTagname(thisStar,"metalYieldMass")
          if (getLength(propertyList) == 1) then
             thisDatum => item(propertyList,0)
             metalYieldCount=metalYieldCount+1
             call extractDataContent(thisDatum,metalYield(metalYieldCount))
             metalYieldMass       (metalYieldCount)=initialMass
             metalYieldMetallicity(metalYieldCount)=metallicity
          end if

          ! Process element yields.
          starHasElements=.false.
          do iElement=1,size(elementYieldCount)
             propertyList => getElementsByTagname(thisStar,"elementYieldMass"//trim(Atomic_Short_Label(iElement)))
             if (getLength(propertyList) == 1) then
                starHasElements=.true.
                thisDatum => item(propertyList,0)
                elementYieldCount(iElement)=elementYieldCount(iElement)+1
                call extractDataContent(thisDatum,elementYield(elementYieldCount(iElement),atomIndexMap(iElement)))
                elementYieldMass       (elementYieldCount(iElement),atomIndexMap(iElement))=initialMass
                elementYieldMetallicity(elementYieldCount(iElement),atomIndexMap(iElement))=metallicity
             end if
          end do

          ! Set any elements that were not included for this star to zero.
          if (starHasElements) then
             do iElement=1,size(elementYieldCount)
                if (elementYieldCount(iElement) < maxval(elementYieldCount) .and. atomIndexMap(iElement) > 0) then
                   elementYieldCount(iElement)=elementYieldCount(iElement)+1
                   elementYieldMass       (elementYieldCount(iElement),atomIndexMap(iElement))=initialMass
                   elementYieldMetallicity(elementYieldCount(iElement),atomIndexMap(iElement))=metallicity
                   elementYield           (elementYieldCount(iElement),atomIndexMap(iElement))=0.0d0
                end if
             end do
          end if

       end do

       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)

    end if

    return
  end subroutine Stellar_Astrophysics_File_Initialize

  double precision function Star_Initial_Mass_File(lifetime,metallicity)
    !% Return the initial mass of a star of given {\tt lifetime} and {\tt metallicity}.
    implicit none
    double precision                         , intent(in   ) :: lifetime                     , metallicity
    type            (interp2dIrregularObject), save          :: interpolationWorkspace
    logical                                  , save          :: resetInterpolation    =.true.
    !$omp threadprivate(interpolationWorkspace,resetInterpolation)

    resetInterpolation=.true.
    Star_Initial_Mass_File=Interpolate_2D_Irregular(stellarLifetime,stellarLifetimeMetallicity,stellarLifetimeMass,lifetime&
         &,metallicity,interpolationWorkspace,reset=resetInterpolation,numberComputePoints=3)
    return
  end function Star_Initial_Mass_File

  double precision function Star_Lifetime_File(initialMass,metallicity)
    !% Return the lifetime of a star (in Gyr) given an {\tt initialMass} and {\tt metallicity}.
    implicit none
    double precision                         , intent(in   ) :: initialMass                  , metallicity
    type            (interp2dIrregularObject), save          :: interpolationWorkspace
    logical                                  , save          :: resetInterpolation    =.true.
    !$omp threadprivate(interpolationWorkspace,resetInterpolation)

    resetInterpolation=.true.
    Star_Lifetime_File=Interpolate_2D_Irregular(stellarLifetimeMass,stellarLifetimeMetallicity,stellarLifetime,initialMass&
         &,metallicity ,interpolationWorkspace,reset=resetInterpolation)

    return
  end function Star_Lifetime_File

  double precision function Star_Ejected_Mass_File(initialMass,metallicity)
    !% Return the mass ejected during the lifetime of a star of given {\tt initialMass} and {\tt metallicity}.
    implicit none
    double precision                         , intent(in   ) :: initialMass                  , metallicity
    type            (interp2dIrregularObject), save          :: interpolationWorkspace
    logical                                  , save          :: resetInterpolation    =.true.
    !$omp threadprivate(interpolationWorkspace,resetInterpolation)

    ! Compute the ejected mass.
    resetInterpolation=.true.
    Star_Ejected_Mass_File=max(Interpolate_2D_Irregular(ejectedMassMass,ejectedMassMetallicity,ejectedMass,initialMass,metallicity&
         &,interpolationWorkspace,reset=resetInterpolation),0.0d0)

    return
  end function Star_Ejected_Mass_File

  double precision function Star_Metal_Yield_Mass_File(initialMass,metallicity,atomIndex)
    !% Return the mass of metals yielded by a star of given {\tt initialMass} and {\tt metallicity}.
    use Memory_Management
    implicit none
    double precision                         , intent(in   )           :: initialMass                  , metallicity
    integer                                  , intent(in   ), optional :: atomIndex
    type            (interp2dIrregularObject), save                    :: interpolationWorkspace
    logical                                  , save                    :: resetInterpolation    =.true.
    !$omp threadprivate(interpolationWorkspace,resetInterpolation)
    integer                                                            :: elementIndex

    if (present(atomIndex)) then
       ! Compute the element mass yield.
       elementIndex=atomIndexMap(atomIndex)
       resetInterpolation=.true.
       Star_Metal_Yield_Mass_File=max(Interpolate_2D_Irregular(elementYieldMass(1:elementYieldCount(atomIndex),elementIndex)&
            &,elementYieldMetallicity(1:elementYieldCount(atomIndex),elementIndex),elementYield(1:elementYieldCount(atomIndex) &
            &,elementIndex),initialMass ,metallicity,interpolationWorkspace,reset =resetInterpolation),0.0d0)
    else
       ! Compute the metal mass yield.
       resetInterpolation=.true.
       Star_Metal_Yield_Mass_File=max(Interpolate_2D_Irregular(metalYieldMass,metalYieldMetallicity,metalYield,initialMass&
            &,metallicity,interpolationWorkspace,reset=resetInterpolation),0.0d0)
    end if
    return
  end function Star_Metal_Yield_Mass_File

end module Stellar_Astrophysics_File
