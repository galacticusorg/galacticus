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

  !% Implements a stellar astrophysics class in which the stellar properties are read from file and interpolated.

  use :: Numerical_Interpolation_2D_Irregular, only : interp2dIrregularObject

  !# <stellarAstrophysics name="stellarAstrophysicsFile">
  !#  <description>
  !#   A stellar astrophysics class which reads properties of individual stars of different initial mass and metallicity from an
  !#   XML file and interpolates in them. The stars can be irregularly spaced in the plane of initial mass and metallicity. The
  !#   XML file should have the following structure:
  !#   \begin{verbatim}
  !#    &lt;stars&gt;
  !#     &lt;star&gt;
  !#       &lt;initialMass&gt;0.6&lt;/initialMass&gt;
  !#       &lt;lifetime&gt;28.19&lt;/lifetime&gt;
  !#       &lt;metallicity&gt;0.0000&lt;/metallicity&gt;
  !#       &lt;ejectedMass&gt;7.65&lt;/ejectedMass&gt;
  !#       &lt;metalYieldMass&gt;0.44435954&lt;/metalYieldMass&gt;
  !#       &lt;elementYieldMassFe&gt;2.2017e-13&lt;/elementYieldMassFe&gt;
  !#       &lt;source&gt;Table 2 of Tumlinson, Shull &amp; Venkatesan (2003, ApJ, 584, 608)&lt;/source&gt;
  !#       &lt;url&gt;http://adsabs.harvard.edu/abs/2003ApJ...584..608T&lt;/url&gt;
  !#     &lt;/star&gt;
  !#     &lt;star&gt;
  !#       .
  !#       .
  !#       .
  !#     &lt;/star&gt;
  !#     .
  !#     .
  !#     .
  !#    &lt;/stars&gt;
  !#   \end{verbatim}
  !#   Each {\normalfont \ttfamily star} element must contain the {\normalfont \ttfamily initialMass} (given in $M_\odot$) and
  !#   {\normalfont \ttfamily metallicity} tags. Other tags are optional. {\normalfont \ttfamily lifetime} gives the lifetime of
  !#   such a star (in Gyr), {\normalfont \ttfamily ejectedMass} gives the total mass (in $M_\odot$) ejected by such a star during
  !#   its lifetime, {\normalfont \ttfamily metalYieldMass} gives the total mass of metals yielded by the star during its lifetime
  !#   while {\normalfont \ttfamily elementYieldMassX} gives the mass of element {\normalfont \ttfamily X} yielded by the star
  !#   during its lifetime. The {\normalfont \ttfamily source} and {\normalfont \ttfamily url} tags are not used, but are strongly
  !#   recommended to provide a reference to the origin of the stellar data.
  !#  </description>
  !# </stellarAstrophysics>
  type, extends(stellarAstrophysicsClass) :: stellarAstrophysicsFile
     !% A stellar astrophysics class in which the stellar properties are read from file and interpolated.
     private
     type            (varying_string         )                              :: fileName
     double precision                         , allocatable, dimension(:  ) :: lifetimeLifetime                 , lifetimeMass                   , &
          &                                                                    lifetimeMetallicity
     double precision                         , allocatable, dimension(:  ) :: massEjectedMassEjected           , massEjectedMass                , &
          &                                                                    massEjectedMetallicity
     double precision                         , allocatable, dimension(:  ) :: yieldMetals                      , yieldMetalsMass                , &
          &                                                                    yieldMetalsMetallicity
     double precision                         , allocatable, dimension(:,:) :: yieldElement                     , yieldElementMass               , &
          &                                                                    yieldElementMetallicity
     integer                                  , allocatable, dimension(:  ) :: atomIndexMap                     , countYieldElement
     integer                                                                :: countElement
     type            (interp2dIrregularObject)                              :: interpolationWorkspaceMassInitial, interpolationWorkspaceLifetime , &
          &                                                                    interpolationWorkspaceMassEjected, interpolationWorkspaceMassYield
     logical                                                                :: interpolationResetMassInitial    , interpolationResetLifetime     , &
          &                                                                    interpolationResetMassEjected    , interpolationResetMassYield    , &
          &                                                                    readDone
  contains
     !# <methods>
     !#   <method description="Read stellar astrophysics data from file." method="read" />
     !# </methods>
     procedure :: massInitial => fileMassInitial
     procedure :: massEjected => fileMassEjected
     procedure :: massYield   => fileMassYield
     procedure :: lifetime    => fileLifetime
     procedure :: read        => fileRead
  end type stellarAstrophysicsFile

  interface stellarAstrophysicsFile
     !% Constructors for the {\normalfont \ttfamily file} stellar astrophysics class.
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface stellarAstrophysicsFile

  ! Current file format version.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function fileConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily file} stellar astrophysics class which takes a parameter list as input.
    use :: Galacticus_Paths, only : galacticusPath, pathTypeDataStatic
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(stellarAstrophysicsFile)                :: self
    type(inputParameters        ), intent(inout) :: parameters
    type(varying_string         )                :: fileName

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <defaultValue>galacticusPath(pathTypeDataStatic)//'stellarAstrophysics/Stellar_Properties_Compilation.xml'</defaultValue>
    !#   <description>The name of the XML file from which to read stellar properties (ejected masses, yields, etc.).</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=stellarAstrophysicsFile(char(fileName))
    !# <inputParametersValidate source="parameters"/>
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName) result(self)
    !% Internal constructor for the {\normalfont \ttfamily file} stellar astrophysics class.
    use :: Atomic_Data      , only : Atomic_Data_Atoms_Count
    use :: Memory_Management, only : allocateArray
    implicit none
    type     (stellarAstrophysicsFile)                :: self
    character(len=*                  ), intent(in   ) :: fileName
    !# <constructorAssign variables="fileName"/>

    ! Allocate array to store number of entries in file for yield of each element.
    call allocateArray(self%countYieldElement,[Atomic_Data_Atoms_Count()])
    call allocateArray(self%atomIndexMap     ,[Atomic_Data_Atoms_Count()])
    self%interpolationResetMassInitial=.true.
    self%interpolationResetLifetime   =.true.
    self%interpolationResetMassEjected=.true.
    self%interpolationResetMassYield  =.true.
    self%readDone                     =.false.
    return
  end function fileConstructorInternal

  subroutine fileRead(self)
    !% Read stellar astrophysics data. This is not done during object construction since it can be slow---we only perform the read if the data is actually needed.
    use :: Atomic_Data      , only : Atomic_Short_Label
    use :: FoX_DOM          , only : destroy                          , node
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: IO_XML           , only : XML_Get_First_Element_By_Tag_Name, XML_Get_Elements_By_Tag_Name, xmlNodeList, extractDataContent => extractDataContentTS, &
         &                           XML_Parse
    use :: Memory_Management, only : allocateArray
    implicit none
    class           (stellarAstrophysicsFile), intent(inout)               :: self
    type            (node                   ), pointer                     :: doc              , datum                   , &
         &                                                                    star
    type            (xmlNodeList            ), allocatable  , dimension(:) :: propertyList     , starList
    integer                                                                :: countMassEjected , countYieldElementMaximum, &
         &                                                                    fileFormatVersion, iElement                , &
         &                                                                    iStar            , ioErr                   , &
         &                                                                    countLifetime    , mapToIndex              , &
         &                                                                    countYieldMetals
    double precision                                                       :: massInitial      , metallicity
    logical                                                                :: starHasElements

    if (self%readDone) return
    !$omp critical (FoX_DOM_Access)
    ! Open the XML file containing stellar properties.
    doc => XML_Parse(char(self%fileName),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Unable to parse stellar properties file'//{introspection:location})
    ! Check the file format version of the file.
    datum => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
    call extractDataContent(datum,fileFormatVersion)
    if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('file format version is out of date'//{introspection:location})
    ! Get a list of all stars.
    call XML_Get_Elements_By_Tag_Name(doc,"star",starList)
    ! Count up number of stars with given properties.
    countLifetime         =0
    countMassEjected      =0
    countYieldMetals      =0
    self%countYieldElement=0
    do iStar=0,size(starList)-1
       star         => starList(iStar)%element
       call XML_Get_Elements_By_Tag_Name(star,"initialMass",propertyList)
       if (size(propertyList) /= 1) call Galacticus_Error_Report('star must have precisely one initial mass'//{introspection:location})
       call XML_Get_Elements_By_Tag_Name(star,"metallicity",propertyList)
       if (size(propertyList) /= 1) call Galacticus_Error_Report('star must have precisely one metallicity' //{introspection:location})
       call XML_Get_Elements_By_Tag_Name(star,"lifetime",propertyList)
       if (size(propertyList) == 1) countLifetime=countLifetime      +1
       if (size(propertyList) >  1) call Galacticus_Error_Report('star has multiple lifetimes'              //{introspection:location})
       call XML_Get_Elements_By_Tag_Name(star,"ejectedMass",propertyList)
       if (size(propertyList) == 1) countMassEjected=countMassEjected+1
       if (size(propertyList) >  1) call Galacticus_Error_Report('star has multiple ejected masses'         //{introspection:location})
       call XML_Get_Elements_By_Tag_Name(star,"metalYieldMass",propertyList)
       if (size(propertyList) == 1) countYieldMetals=countYieldMetals+1
       if (size(propertyList) >  1) call Galacticus_Error_Report('star has multiple metal yield masses'     //{introspection:location})
       do iElement=1,size(self%countYieldElement)
          call XML_Get_Elements_By_Tag_Name(star,"elementYieldMass"//trim(Atomic_Short_Label(iElement)),propertyList)
          if (size(propertyList) == 1) self%countYieldElement(iElement)=self%countYieldElement(iElement)+1
          if (size(propertyList) >  1) call Galacticus_Error_Report('star has multiple element yield masses'//{introspection:location})
       end do
    end do
    ! Find number of elements for which some yield data is available.
    self%countElement       =count (self%countYieldElement > 0)
    countYieldElementMaximum=maxval(self%countYieldElement    )
    ! Create mapping of atomic index to our array space.
    mapToIndex=0
    do iElement=1,size(self%countYieldElement)
       if (self%countYieldElement(iElement) > 0) then
          mapToIndex=mapToIndex+1
          self%atomIndexMap(iElement)=mapToIndex
       else
          self%atomIndexMap(iElement)=-1
       end if
    end do
    ! Allocate arrays to store stellar properties.
    call allocateArray(self%lifetimeLifetime          ,[countLifetime                             ])
    call allocateArray(self%lifetimeMass              ,[countLifetime                             ])
    call allocateArray(self%lifetimeMetallicity       ,[countLifetime                             ])
    call allocateArray(self%massEjectedMassEjected    ,[countMassEjected                          ])
    call allocateArray(self%massEjectedMass           ,[countMassEjected                          ])
    call allocateArray(self%massEjectedMetallicity    ,[countMassEjected                          ])
    call allocateArray(self%yieldMetals               ,[countYieldMetals                          ])
    call allocateArray(self%yieldMetalsMass           ,[countYieldMetals                          ])
    call allocateArray(self%yieldMetalsMetallicity    ,[countYieldMetals                          ])
    call allocateArray(self%yieldElement              ,[countYieldElementMaximum,self%countElement])
    call allocateArray(self%yieldElementMass          ,[countYieldElementMaximum,self%countElement])
    call allocateArray(self%yieldElementMetallicity   ,[countYieldElementMaximum,self%countElement])
    ! Loop over stars to process their properties.
    countLifetime         =0
    countMassEjected      =0
    countYieldMetals      =0
    self%countYieldElement=0
    do iStar=0,size(starList)-1
       star         => starList(iStar)%element
       call XML_Get_Elements_By_Tag_Name(star,"initialMass",propertyList)
       datum        => propertyList(0)%element
       call extractDataContent(datum,massInitial)
       call XML_Get_Elements_By_Tag_Name(star,"metallicity",propertyList)
       datum        => propertyList(0)%element
       call extractDataContent(datum,metallicity)
       ! Process stellar lifetimes.
       call XML_Get_Elements_By_Tag_Name(star,"lifetime",propertyList)
       if (size(propertyList) == 1) then
          datum => propertyList(0)%element
          countLifetime=countLifetime+1
          call extractDataContent(datum,self%lifetimeLifetime      (countLifetime   ))
          self%lifetimeMass       (countLifetime)=massInitial
          self%lifetimeMetallicity(countLifetime)=metallicity
       end if
       ! Process ejected masses.
       call XML_Get_Elements_By_Tag_Name(star,"ejectedMass",propertyList)
       if (size(propertyList) == 1) then
          datum => propertyList(0)%element
          countMassEjected=countMassEjected+1
          call extractDataContent(datum,self%massEjectedMassEjected(countMassEjected))
          self%massEjectedMass       (countMassEjected)=massInitial
          self%massEjectedMetallicity(countMassEjected)=metallicity
       end if
       ! Process metal yields.
       call XML_Get_Elements_By_Tag_Name(star,"metalYieldMass",propertyList)
       if (size(propertyList) == 1) then
          datum => propertyList(0)%element
          countYieldMetals=countYieldMetals+1
          call extractDataContent(datum,self%yieldMetals           (countYieldMetals))
          self%yieldMetalsMass       (countYieldMetals)=massInitial
          self%yieldMetalsMetallicity(countYieldMetals)=metallicity
       end if
       ! Process element yields.
       starHasElements=.false.
       do iElement=1,size(self%countYieldElement)
          call XML_Get_Elements_By_Tag_Name(star,"elementYieldMass"//trim(Atomic_Short_Label(iElement)),propertyList)
          if (size(propertyList) == 1) then
             starHasElements=.true.
             datum => propertyList(0)%element
             self%countYieldElement(iElement)=self%countYieldElement(iElement)+1
             call extractDataContent(datum,self%yieldElement(self%countYieldElement(iElement),self%atomIndexMap(iElement)))
             self%yieldElementMass       (self%countYieldElement(iElement),self%atomIndexMap(iElement))=massInitial
             self%yieldElementMetallicity(self%countYieldElement(iElement),self%atomIndexMap(iElement))=metallicity
          end if
       end do
       ! Set any elements that were not included for this star to zero.
       if (starHasElements) then
          do iElement=1,size(self%countYieldElement)
             if (self%countYieldElement(iElement) < maxval(self%countYieldElement) .and. self%atomIndexMap(iElement) > 0) then
                self%countYieldElement(iElement)=self%countYieldElement(iElement)+1
                self%yieldElementMass       (self%countYieldElement(iElement),self%atomIndexMap(iElement))=massInitial
                self%yieldElementMetallicity(self%countYieldElement(iElement),self%atomIndexMap(iElement))=metallicity
                self%yieldElement           (self%countYieldElement(iElement),self%atomIndexMap(iElement))=0.0d0
             end if
          end do
       end if
    end do
    ! Destroy the document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    self%readDone=.true.
    return
  end subroutine fileRead

  double precision function fileMassInitial(self,lifetime,metallicity)
    !% Return the initial mass of a star of given {\normalfont \ttfamily lifetime} and {\normalfont \ttfamily metallicity}.
    use :: Numerical_Interpolation_2D_Irregular, only : Interpolate_2D_Irregular
    implicit none
    class           (stellarAstrophysicsFile), intent(inout) :: self
    double precision                         , intent(in   ) :: lifetime, metallicity

    call self%read()
    fileMassInitial=Interpolate_2D_Irregular(                                                            &
         &                                                       self%lifetimeLifetime                 , &
         &                                                       self%lifetimeMetallicity              , &
         &                                                       self%lifetimeMass                     , &
         &                                                            lifetime                         , &
         &                                                            metallicity                      , &
         &                                                       self%interpolationWorkspaceMassInitial, &
         &                                   reset              =self%interpolationResetMassInitial    , &
         &                                   numberComputePoints=     3                                  &
         &                                  )
    return
  end function fileMassInitial

  double precision function fileLifetime(self,massInitial,metallicity)
    !% Return the lifetime of a star (in Gyr) given an {\normalfont \ttfamily massInitial} and {\normalfont \ttfamily metallicity}.
    use :: Numerical_Interpolation_2D_Irregular, only : Interpolate_2D_Irregular
    implicit none
    class           (stellarAstrophysicsFile), intent(inout) :: self
    double precision                         , intent(in   ) :: massInitial, metallicity

    call self%read()
    fileLifetime=Interpolate_2D_Irregular(                                           &
         &                                      self%lifetimeMass                  , &
         &                                      self%lifetimeMetallicity           , &
         &                                      self%lifetimeLifetime              , &
         &                                           massInitial                   , &
         &                                           metallicity                   , &
         &                                      self%interpolationWorkspaceLifetime, &
         &                                reset=self%interpolationResetLifetime      &
         &                               )
    return
  end function fileLifetime

  double precision function fileMassEjected(self,massInitial,metallicity)
    !% Return the mass ejected during the lifetime of a star of given {\normalfont \ttfamily massInitial} and {\normalfont \ttfamily metallicity}.
    use :: Numerical_Interpolation_2D_Irregular, only : Interpolate_2D_Irregular
    implicit none
    class           (stellarAstrophysicsFile), intent(inout) :: self
    double precision                         , intent(in   ) :: massInitial, metallicity

    call self%read()
    fileMassEjected=max(                                                                       &
         &              Interpolate_2D_Irregular(                                              &
         &                                             self%massEjectedMass                  , &
         &                                             self%massEjectedMetallicity           , &
         &                                             self%massEjectedMassEjected           , &
         &                                                  massInitial                      , &
         &                                                  metallicity                      , &
         &                                             self%interpolationWorkspaceMassEjected, &
         &                                       reset=self%interpolationResetMassEjected      &
         &                                      )                                            , &
         &               0.0d0                                                                 &
         &             )
    return
  end function fileMassEjected

  double precision function fileMassYield(self,massInitial,metallicity,atomIndex)
    !% Return the mass of metals yielded by a star of given {\normalfont \ttfamily massInitial} and {\normalfont \ttfamily metallicity}.
    use :: Numerical_Interpolation_2D_Irregular, only : Interpolate_2D_Irregular
    implicit none
    class           (stellarAstrophysicsFile), intent(inout)           :: self
    double precision                         , intent(in   )           :: massInitial , metallicity
    integer                                  , intent(in   ), optional :: atomIndex
    integer                                                            :: elementIndex

    call self%read()
    if (present(atomIndex)) then
       ! Compute the element mass yield.
       elementIndex =self%atomIndexMap(atomIndex)
       fileMassYield=max(                                                                                                                       &
            &            Interpolate_2D_Irregular(                                                                                              &
            &                                           self%yieldElementMass               (1:self%countYieldElement(atomIndex),elementIndex), &
            &                                           self%yieldElementMetallicity        (1:self%countYieldElement(atomIndex),elementIndex), &
            &                                           self%yieldElement                   (1:self%countYieldElement(atomIndex),elementIndex), &
            &                                                massInitial                                                                      , &
            &                                                metallicity                                                                      , &
            &                                           self%interpolationWorkspaceMassYield                                                  , &
            &                                     reset=self%interpolationResetMassYield                                                        &
            &                                    )                                                                                            , &
            &            0.0d0                                                                                                                  &
            &           )
    else
       ! Compute the metal mass yield.
       fileMassYield=max(                                                                     &
            &            Interpolate_2D_Irregular(                                            &
            &                                           self%yieldMetalsMass                , &
            &                                           self%yieldMetalsMetallicity         , &
            &                                           self%yieldMetals                    , &
            &                                                massInitial                    , &
            &                                                metallicity                    , &
            &                                           self%interpolationWorkspaceMassYield, &
            &                                     reset=self%interpolationResetMassYield      &
            &                                    )                                          , &
            &             0.0d0                                                               &
            &            )
    end if
    return
  end function fileMassYield
