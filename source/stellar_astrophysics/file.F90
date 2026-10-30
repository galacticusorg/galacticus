!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a stellar astrophysics class in which the stellar properties are read from file and interpolated.
  !!}

  use :: Numerical_Interpolation_2D_Irregular, only : interpolator2DIrregular

  type :: stellarAstrophysicsFileTable
     !!{RST
     A stellar property tabulated by metallicity. At each tabulated metallicity, ``mass(1:countMass(i),i)`` holds the initial
     masses of the stars (in ascending order) and ``value(1:countMass(i),i)`` the corresponding values of the property.
     !!}
     integer                                       :: countMetallicity
     double precision, allocatable, dimension(:  ) :: metallicity
     integer         , allocatable, dimension(:  ) :: countMass
     double precision, allocatable, dimension(:,:) :: mass            , value
  end type stellarAstrophysicsFileTable

  !![
  <stellarAstrophysics name="stellarAstrophysicsFile" docformat="rst">
   <description>
   A stellar astrophysics class which reads properties of individual stars of different initial mass and metallicity from an XML file and interpolates in them. The stars can be irregularly spaced in the plane of initial mass and metallicity. The XML file should have the following structure:

   .. code-block:: none

       &lt;stars&gt;
        &lt;star&gt;
          &lt;initialMass&gt;0.6&lt;/initialMass&gt;
          &lt;lifetime&gt;28.19&lt;/lifetime&gt;
          &lt;metallicity&gt;0.0000&lt;/metallicity&gt;
          &lt;ejectedMass&gt;7.65&lt;/ejectedMass&gt;
          &lt;metalYieldMass&gt;0.44435954&lt;/metalYieldMass&gt;
          &lt;elementYieldMassFe&gt;2.2017e-13&lt;/elementYieldMassFe&gt;
          &lt;source&gt;Table 2 of Tumlinson, Shull &amp; Venkatesan (2003, ApJ, 584, 608)&lt;/source&gt;
          &lt;url&gt;http://adsabs.harvard.edu/abs/2003ApJ...584..608T&lt;/url&gt;
        &lt;/star&gt;
        &lt;star&gt;
          .
          .
          .
        &lt;/star&gt;
        .
        .
        .
       &lt;/stars&gt;

   Each ``star`` element must contain the ``initialMass`` (given in :math:`\mathrm{M}_\odot`) and ``metallicity`` tags. Other tags are optional. ``lifetime`` gives the lifetime of such a star (in Gyr), ``ejectedMass`` gives the total mass (in :math:`\mathrm{M}_\odot`) ejected by such a star during its lifetime, ``metalYieldMass`` gives the total mass of metals yielded by the star during its lifetime while ``elementYieldMassX`` gives the mass of element ``X`` yielded by the star during its lifetime. The ``source`` and ``url`` tags are not used, but are strongly recommended to provide a reference to the origin of the stellar data.

   The minimum initial mass of a star which ends its life as a Type II supernova is set by the ``[initialMassForSupernovaeTypeII]`` parameter.

   A star which provides a yield for any element is assumed to yield zero mass of every element (for which any star provides a yield) that it does not list.

   Lifetimes (and their inverse, the initial mass of a star of given lifetime) are interpolated in the plane of initial mass and metallicity using a bivariate interpolation on the irregularly spaced points. Ejected masses, metal yields, and element yields are typically compiled from different sources which each cover a different region of that plane, such that the points do not cover it uniformly---for example, the asymptotic giant branch and massive star yields in the standard compilation are tabulated at different sets of metallicities. A bivariate interpolation across such gaps can return arbitrary values. These properties are therefore interpolated linearly in initial mass between the stars tabulated at each metallicity, and then linearly in metallicity between the two tabulated metallicities nearest to (and bracketing) the requested metallicity, considering only those metallicities at which the tabulated stars extend down to the requested initial mass. If the requested metallicity lies outside the range of such metallicities, the value at the nearest is used. If no tabulated metallicity extends down to the requested initial mass, the property is zero. Above the most massive star tabulated at a given metallicity, the property at that metallicity is taken to be that of the most massive star scaled in proportion to initial mass (i.e. its value per unit initial mass is held fixed). This keeps each metallicity's own models in use where its tabulation stops short of others---for example, the Limongi &amp; Chieffi (2018) models extend to 120 :math:`\mathrm{M}_\odot` at their highest metallicity but only to 80 :math:`\mathrm{M}_\odot` at lower metallicities---rather than taking the property from a tabulation at a very different metallicity. It also means that yields tabulated only at zero metallicity (e.g. those of pair instability supernovae) do not apply at higher metallicities. No such extension is made below the least massive star tabulated at a metallicity, since that would extrapolate one class of stellar model into the mass range of another. A metallicity whose tabulation covers a narrower range of masses than those of its neighbors therefore enters the interpolation only at its least massive star, and properties can change abruptly at that mass. (In ``stellarPropertiesCompilationKarakas2010Doherty2014LimongiChieffi2018.xml``, for example, the only models at :math:`Z=0.001` are super-AGB stars of 6.5 :math:`\mathrm{M}_\odot` and above, so yields at nearby metallicities change at that mass.)

   Yields change discontinuously at ``[initialMassForSupernovaeTypeII]``, where the stars that end their lives on the asymptotic giant branch give way to those that end as supernovae, and the two populations are usually tabulated from different sources with a gap between them (from 5 to 9 :math:`\mathrm{M}_\odot` in the standard compilation). Ejected masses and yields are therefore never interpolated across that mass: at each metallicity only the tabulated stars on the same side of it as the requested initial mass are used. Below ``[initialMassForSupernovaeTypeII]``, the tabulation at each metallicity is extended upward to it as described above. A star lying between ``[initialMassForSupernovaeTypeII]`` and the least massive tabulated star above it (at any metallicity) is assigned the ejected mass and yields per unit initial mass of a star of that least massive tabulated mass.
   </description>
   <runTimeFileDependencies paths="fileName"/>
  </stellarAstrophysics>
  !!]
  type, extends(stellarAstrophysicsClass) :: stellarAstrophysicsFile
     !!{RST
     A stellar astrophysics class in which the stellar properties are read from file and interpolated.
     !!}
     private
     type            (varying_string              )                              :: fileName
     double precision                              , allocatable, dimension(:  ) :: lifetimeLifetime                 , lifetimeMass                   , &
          &                                                                         lifetimeMetallicity
     type            (stellarAstrophysicsFileTable)                              :: tableMassEjected                 , tableYieldMetals
     ! One table per element for which yields are available.
     type            (stellarAstrophysicsFileTable), allocatable, dimension(:  ) :: tableYieldElement
     integer                                       , allocatable, dimension(:  ) :: atomIndexMap                     , countYieldElement
     integer                                                                     :: countElement
     double precision                                                            :: initialMassForSupernovaeTypeII
     type            (interpolator2DIrregular     )                              :: interpolationWorkspaceMassInitial, interpolationWorkspaceLifetime
     logical                                                                     :: interpolationResetMassInitial    , interpolationResetLifetime    , &
          &                                                                         readDone
  contains
     !![
     <methods docformat="rst">
       <method description="Read stellar astrophysics data from file." method="read" />
     </methods>
     !!]
     procedure :: massInitial                 => fileMassInitial
     procedure :: massEjected                 => fileMassEjected
     procedure :: massYield                   => fileMassYield
     procedure :: lifetime                    => fileLifetime
     procedure :: massInitialSupernovaeTypeII => fileMassInitialSupernovaeTypeII
     procedure :: read                        => fileRead
  end type stellarAstrophysicsFile

  interface stellarAstrophysicsFile
     !!{RST
     Constructors for the :galacticus-class:`stellarAstrophysicsFile` stellar astrophysics class.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface stellarAstrophysicsFile

  ! Current file format version.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function fileConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`stellarAstrophysicsFile` stellar astrophysics class which takes a parameter list as input.
    !!}
    use :: Input_Paths     , only : inputPath     , pathTypeDataStatic
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarAstrophysicsFile)                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    type            (varying_string         )                :: fileName
    double precision                                         :: initialMassForSupernovaeTypeII

    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <defaultValue>inputPath(pathTypeDataStatic)//'stellarAstrophysics/stellarPropertiesCompilationStandard.xml'</defaultValue>
      <description>
      The name of the XML file from which to read stellar properties (ejected masses, yields, etc.).
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>initialMassForSupernovaeTypeII</name>
      <defaultValue>8.0d0</defaultValue>
      <description>
      The minimum initial mass (in :math:`\mathrm{M}_\odot`) of a star which ends its life as a Type II supernova. Ejected masses and yields are not interpolated across this mass.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=stellarAstrophysicsFile(char(fileName),initialMassForSupernovaeTypeII)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName,initialMassForSupernovaeTypeII) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`stellarAstrophysicsFile` stellar astrophysics class.
    !!}
    use :: Atomic_Data      , only : Atomic_Data_Atoms_Count
    implicit none
    type            (stellarAstrophysicsFile)                :: self
    character       (len=*                  ), intent(in   ) :: fileName
    double precision                         , intent(in   ) :: initialMassForSupernovaeTypeII
    !![
    <constructorAssign variables="fileName, initialMassForSupernovaeTypeII"/>
    !!]

    ! Allocate array to store number of entries in file for yield of each element.
    allocate(self%countYieldElement(Atomic_Data_Atoms_Count()))
    allocate(self%atomIndexMap     (Atomic_Data_Atoms_Count()))
    self%interpolationResetMassInitial=.true.
    self%interpolationResetLifetime   =.true.
    self%readDone                     =.false.
    return
  end function fileConstructorInternal

  subroutine fileRead(self)
    !!{RST
    Read stellar astrophysics data. This is not done during object construction since it can be slow---we only perform the read if the data is actually needed.
    !!}
    use :: Atomic_Data   , only : Atomic_Short_Label
    use :: FoX_DOM       , only : destroy                          , node                        , extractDataContent
    use :: Error         , only : Error_Report
    use :: IO_XML        , only : XML_Get_First_Element_By_Tag_Name, XML_Get_Elements_By_Tag_Name, xmlNodeList       ,  XML_Parse
    implicit none
    class           (stellarAstrophysicsFile), intent(inout)                 :: self
    type            (node                   ), pointer                       :: doc                    , datum                  , &
         &                                                                      star
    type            (xmlNodeList            ), allocatable  , dimension(:  ) :: propertyList           , starList
    double precision                         , allocatable  , dimension(:  ) :: massEjectedMassEjected , massEjectedMass        , &
         &                                                                      massEjectedMetallicity , yieldMetals            , &
         &                                                                      yieldMetalsMass        , yieldMetalsMetallicity
    double precision                         , allocatable  , dimension(:,:) :: yieldElement           , yieldElementMass       , &
         &                                                                      yieldElementMetallicity
    integer                                                                  :: countMassEjected       , countYieldElementStars , &
         &                                                                      fileFormatVersion      , iElement               , &
         &                                                                      iStar                  , ioErr                  , &
         &                                                                      countLifetime          , mapToIndex             , &
         &                                                                      countYieldMetals
    double precision                                                         :: massInitial            , metallicity
    logical                                                                  :: starHasElements

    if (self%readDone) return
    !$omp critical (FoX_DOM_Access)
    ! Open the XML file containing stellar properties.
    doc => XML_Parse(self%fileName,iostat=ioErr)
    if (ioErr /= 0) call Error_Report('Unable to parse stellar properties file'//{introspection:location})
    ! Check the file format version of the file.
    datum => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
    call extractDataContent(datum,fileFormatVersion)
    if (fileFormatVersion /= fileFormatVersionCurrent) call Error_Report('file format version is out of date'//{introspection:location})
    ! Get a list of all stars.
    call XML_Get_Elements_By_Tag_Name(doc,"star",starList)
    ! Count up number of stars with given properties.
    countLifetime         =0
    countMassEjected      =0
    countYieldMetals      =0
    countYieldElementStars=0
    self%countYieldElement=0
    do iStar=0,size(starList)-1
       star         => starList(iStar)%element
       call XML_Get_Elements_By_Tag_Name(star,"initialMass",propertyList)
       if (size(propertyList) /= 1) call Error_Report('star must have precisely one initial mass'//{introspection:location})
       call XML_Get_Elements_By_Tag_Name(star,"metallicity",propertyList)
       if (size(propertyList) /= 1) call Error_Report('star must have precisely one metallicity' //{introspection:location})
       call XML_Get_Elements_By_Tag_Name(star,"lifetime",propertyList)
       if (size(propertyList) == 1) countLifetime=countLifetime      +1
       if (size(propertyList) >  1) call Error_Report('star has multiple lifetimes'              //{introspection:location})
       call XML_Get_Elements_By_Tag_Name(star,"ejectedMass",propertyList)
       if (size(propertyList) == 1) countMassEjected=countMassEjected+1
       if (size(propertyList) >  1) call Error_Report('star has multiple ejected masses'         //{introspection:location})
       call XML_Get_Elements_By_Tag_Name(star,"metalYieldMass",propertyList)
       if (size(propertyList) == 1) countYieldMetals=countYieldMetals+1
       if (size(propertyList) >  1) call Error_Report('star has multiple metal yield masses'     //{introspection:location})
       starHasElements=.false.
       do iElement=1,size(self%countYieldElement)
          call XML_Get_Elements_By_Tag_Name(star,"elementYieldMass"//trim(Atomic_Short_Label(iElement)),propertyList)
          if (size(propertyList) == 1) then
             starHasElements                 =.true.
             self%countYieldElement(iElement)=self%countYieldElement(iElement)+1
          end if
          if (size(propertyList) >  1) call Error_Report('star has multiple element yield masses'//{introspection:location})
       end do
       if (starHasElements) countYieldElementStars=countYieldElementStars+1
    end do
    ! Find number of elements for which some yield data is available.
    self%countElement=count(self%countYieldElement > 0)
    ! Validate.
    if (countLifetime    <= 0) call Error_Report('star compilation provides no lifetimes'     //{introspection:location})
    if (countMassEjected <= 0) call Error_Report('star compilation provides no ejected masses'//{introspection:location})
    if (countYieldMetals <= 0) call Error_Report('star compilation provides no metal yields'  //{introspection:location})
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
    ! Allocate arrays to store stellar properties. Every star which provides any element yield contributes a point to the
    ! table of every element (with zero yield for those elements that it does not list), so the element tables are sized by
    ! the number of such stars.
    allocate(self%lifetimeLifetime       (countLifetime                           ))
    allocate(self%lifetimeMass           (countLifetime                           ))
    allocate(self%lifetimeMetallicity    (countLifetime                           ))
    allocate(     massEjectedMassEjected (countMassEjected                        ))
    allocate(     massEjectedMass        (countMassEjected                        ))
    allocate(     massEjectedMetallicity (countMassEjected                        ))
    allocate(     yieldMetals            (countYieldMetals                        ))
    allocate(     yieldMetalsMass        (countYieldMetals                        ))
    allocate(     yieldMetalsMetallicity (countYieldMetals                        ))
    allocate(     yieldElement           (countYieldElementStars,self%countElement))
    allocate(     yieldElementMass       (countYieldElementStars,self%countElement))
    allocate(     yieldElementMetallicity(countYieldElementStars,self%countElement))
    ! Loop over stars to process their properties.
    countLifetime         =0
    countMassEjected      =0
    countYieldMetals      =0
    countYieldElementStars=0
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
          call extractDataContent(datum,massEjectedMassEjected(countMassEjected))
          massEjectedMass       (countMassEjected)=massInitial
          massEjectedMetallicity(countMassEjected)=metallicity
       end if
       ! Process metal yields.
       call XML_Get_Elements_By_Tag_Name(star,"metalYieldMass",propertyList)
       if (size(propertyList) == 1) then
          datum => propertyList(0)%element
          countYieldMetals=countYieldMetals+1
          call extractDataContent(datum,yieldMetals(countYieldMetals))
          yieldMetalsMass       (countYieldMetals)=massInitial
          yieldMetalsMetallicity(countYieldMetals)=metallicity
       end if
       ! Process element yields. Elements not listed for a star which lists any element are assigned zero yield.
       starHasElements=.false.
       do iElement=1,size(self%countYieldElement)
          call XML_Get_Elements_By_Tag_Name(star,"elementYieldMass"//trim(Atomic_Short_Label(iElement)),propertyList)
          if (size(propertyList) == 1) then
             if (.not.starHasElements) then
                starHasElements       =.true.
                countYieldElementStars=countYieldElementStars+1
                yieldElement           (countYieldElementStars,:)=0.0d0
                yieldElementMass       (countYieldElementStars,:)=massInitial
                yieldElementMetallicity(countYieldElementStars,:)=metallicity
             end if
             datum => propertyList(0)%element
             call extractDataContent(datum,yieldElement(countYieldElementStars,self%atomIndexMap(iElement)))
          end if
       end do
    end do
    ! Destroy the document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    ! Check that the lifetimes span a range of metallicity. Lifetimes tabulated at a single metallicity cannot be interpolated
    ! in the (mass, metallicity) plane, and requesting them at any other metallicity silently returns nonsense rather than
    ! failing, so it is caught here instead.
    call fileValidateMetallicities(self%fileName,'stellar lifetimes',self%lifetimeMetallicity)
    ! Construct tables of ejected masses and yields.
    call fileTableBuild(self%fileName,'ejected masses'       ,self%tableMassEjected,massEjectedMass,massEjectedMetallicity,massEjectedMassEjected)
    call fileTableBuild(self%fileName,'the total metal yield',self%tableYieldMetals,yieldMetalsMass,yieldMetalsMetallicity,yieldMetals           )
    allocate(self%tableYieldElement(self%countElement))
    do iElement=1,size(self%countYieldElement)
       if (self%atomIndexMap(iElement) > 0)                                                                           &
            & call fileTableBuild(                                                                                    &
            &                     self%fileName                                                                     , &
            &                     'the yield of '//trim(Atomic_Short_Label(iElement))                               , &
            &                     self%tableYieldElement      (                         self%atomIndexMap(iElement)), &
            &                          yieldElementMass       (1:countYieldElementStars,self%atomIndexMap(iElement)), &
            &                          yieldElementMetallicity(1:countYieldElementStars,self%atomIndexMap(iElement)), &
            &                          yieldElement           (1:countYieldElementStars,self%atomIndexMap(iElement))  &
            &                    )
    end do
    self%readDone=.true.
    return
  end subroutine fileRead

  subroutine fileTableBuild(fileName,label,table,mass,metallicity,value)
    !!{RST
    Organize the stellar property ``value``, tabulated for stars of given initial ``mass`` and ``metallicity``, into a
    ``table`` ordered by metallicity, and by initial mass at each metallicity.
    !!}
    use            :: Error             , only : Error_Report
    use            :: ISO_Varying_String, only : char
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: Sorting           , only : sortIndex
    implicit none
    type            (varying_string              ), intent(in   )               :: fileName
    character       (len=*                       ), intent(in   )               :: label
    type            (stellarAstrophysicsFileTable), intent(  out)               :: table
    double precision                              , intent(in   ), dimension(:) :: mass             , metallicity, &
         &                                                                         value
    integer         (c_size_t                    ), allocatable  , dimension(:) :: order
    double precision                              , allocatable  , dimension(:) :: metallicityUnique
    integer                                                                     :: i                , j
    character       (len=32                      )                              :: label_

    ! Find the unique metallicities, in ascending order.
    order=sortIndex(metallicity)
    allocate(metallicityUnique(size(metallicity)))
    table%countMetallicity=0
    do i=1,size(order)
       if (table%countMetallicity > 0) then
          if (metallicity(order(i)) == metallicityUnique(table%countMetallicity)) cycle
       end if
       table%countMetallicity                   =table%countMetallicity+1
       metallicityUnique(table%countMetallicity)=metallicity(order(i))
    end do
    table%metallicity=metallicityUnique(1:table%countMetallicity)
    ! Gather the stars at each metallicity, in ascending order of initial mass.
    allocate(table%countMass(table%countMetallicity                 ))
    allocate(table%mass     (size(mass),table%countMetallicity))
    allocate(table%value    (size(mass),table%countMetallicity))
    table%countMass=0
    order          =sortIndex(mass)
    do i=1,size(order)
       j=findloc(table%metallicity,metallicity(order(i)),dim=1)
       if (table%countMass(j) > 0) then
          if (mass(order(i)) == table%mass(table%countMass(j),j)) then
             write (label_,'(e12.6,", ",e12.6)') mass(order(i)),metallicity(order(i))
             call Error_Report(                                                                                             &
                  &            'multiple values of '//label//' are given for the star of (initial mass, metallicity) = ('// &
                  &            trim(adjustl(label_))//') in stellar properties file "'//char(fileName)//'"'              // &
                  &            {introspection:location}                                                                     &
                  &           )
          end if
       end if
       table%countMass(                  j)=table%countMass(j)+1
       table%mass     (table%countMass(j),j)=mass (order(i))
       table%value    (table%countMass(j),j)=value(order(i))
    end do
    return
  end subroutine fileTableBuild

  double precision function fileTableInterpolate(table,massInitial,metallicity,massThreshold) result(value)
    !!{RST
    Interpolate in a ``table`` of a stellar property to the given ``massInitial`` and ``metallicity``. The property is
    interpolated linearly in initial mass at each tabulated metallicity, and then linearly in metallicity between the nearest
    tabulated metallicities, below and above the requested metallicity, at which the tabulated stars extend down to the
    requested initial mass. If such a metallicity exists on only one side of the requested metallicity, the value there is
    used. If none exists, zero is returned. Above the most massive star tabulated at a metallicity, the property at that
    metallicity is taken to be that of the most massive star, scaled in proportion to initial mass (i.e. the value per unit
    initial mass is held fixed).

    If ``massThreshold`` is given, the stars below it and those at or above it are treated as separate tables, and the
    interpolation in initial mass never crosses it. Stars below the threshold are therefore extended up to it as described
    above. A star lying between the threshold and the least massive tabulated star above it (at any metallicity) is assigned
    the value per unit initial mass of a star of that least massive tabulated mass.
    !!}
    implicit none
    type            (stellarAstrophysicsFileTable), intent(in   )           :: table
    double precision                              , intent(in   )           :: massInitial  , metallicity
    double precision                              , intent(in   ), optional :: massThreshold
    integer                                                                 :: i            , indexLow    , &
         &                                                                     indexHigh    , indexMinimum, &
         &                                                                     indexMaximum
    double precision                                                        :: factor       , massEvaluate

    ! Determine the initial mass at which to evaluate the property.
    massEvaluate=massInitial
    if (present(massThreshold)) then
       if (massInitial >= massThreshold) then
          ! Find the least massive tabulated star at or above the threshold, at any metallicity. Between the threshold and
          ! such a star the property is evaluated at the mass of that star. (Below the threshold no such adjustment is needed,
          ! as the tabulation at each metallicity is extended upward to the threshold.)
          massEvaluate=+huge(0.0d0)
          do i=1,table%countMetallicity
             call fileTableIndexRange(table,i,massInitial,massThreshold,indexMinimum,indexMaximum)
             if (indexMaximum >= indexMinimum) massEvaluate=min(massEvaluate,table%mass(indexMinimum,i))
          end do
          if (massEvaluate == +huge(0.0d0)) then
             value=0.0d0
             return
          end if
          massEvaluate=max(massInitial,massEvaluate)
       end if
    end if
    ! Find the nearest metallicities, below (or equal to) and above the requested metallicity, at which the tabulated stars
    ! (on the same side of any threshold) extend down to the initial mass.
    indexLow =0
    indexHigh=0
    do i=1,table%countMetallicity
       call fileTableIndexRange(table,i,massInitial,massThreshold,indexMinimum,indexMaximum)
       if (indexMaximum < indexMinimum) cycle
       if (massEvaluate < table%mass(indexMinimum,i)) cycle
       if (table%metallicity(i) <= metallicity) then
          indexLow =i
       else
          indexHigh=i
          exit
       end if
    end do
    if      (indexLow == 0 .and. indexHigh == 0) then
       value=0.0d0
    else if (indexLow == 0                     ) then
       value=fileTableInterpolateMass(table,massEvaluate,indexHigh,massInitial,massThreshold)
    else if (                    indexHigh == 0) then
       value=fileTableInterpolateMass(table,massEvaluate,indexLow ,massInitial,massThreshold)
    else
       factor=+(      metallicity           -table%metallicity(indexLow)) &
            & /(table%metallicity(indexHigh)-table%metallicity(indexLow))
       value =+(1.0d0-factor)*fileTableInterpolateMass(table,massEvaluate,indexLow ,massInitial,massThreshold) &
            & +       factor *fileTableInterpolateMass(table,massEvaluate,indexHigh,massInitial,massThreshold)
    end if
    ! Scale to the requested initial mass.
    value  =+value        &
         &  *massInitial  &
         &  /massEvaluate
    return
  end function fileTableInterpolate

  subroutine fileTableIndexRange(table,indexMetallicity,massInitial,massThreshold,indexMinimum,indexMaximum)
    !!{RST
    Find the range of indices, ``indexMinimum`` to ``indexMaximum``, of the stars tabulated in ``table`` at the metallicity
    with index ``indexMetallicity`` which lie on the same side of ``massThreshold`` (if given) as ``massInitial``. The range is
    empty (``indexMaximum`` < ``indexMinimum``) if there are no such stars.
    !!}
    implicit none
    type            (stellarAstrophysicsFileTable), intent(in   )           :: table
    integer                                       , intent(in   )           :: indexMetallicity
    double precision                              , intent(in   )           :: massInitial
    double precision                              , intent(in   ), optional :: massThreshold
    integer                                       , intent(  out)           :: indexMinimum    , indexMaximum
    integer                                                                 :: countBelow

    indexMinimum=1
    indexMaximum=table%countMass(indexMetallicity)
    if (present(massThreshold)) then
       ! Stars are sorted by initial mass, so those below the threshold precede those at or above it.
       countBelow=count(table%mass(1:table%countMass(indexMetallicity),indexMetallicity) < massThreshold)
       if (massInitial < massThreshold) then
          indexMaximum=countBelow
       else
          indexMinimum=countBelow+1
       end if
    end if
    return
  end subroutine fileTableIndexRange

  double precision function fileTableInterpolateMass(table,massEvaluate,indexMetallicity,massInitial,massThreshold) result(value)
    !!{RST
    Interpolate linearly in initial mass in a ``table`` of a stellar property to ``massEvaluate`` at the tabulated metallicity
    with index ``indexMetallicity``, using only those tabulated stars which lie on the same side of ``massThreshold`` (if
    given) as ``massInitial``. ``massEvaluate`` must not lie below the least massive of those stars. Above the most massive
    of them, the value for that star is scaled in proportion to ``massEvaluate``.
    !!}
    implicit none
    type            (stellarAstrophysicsFileTable), intent(in   )           :: table
    double precision                              , intent(in   )           :: massEvaluate    , massInitial
    integer                                       , intent(in   )           :: indexMetallicity
    double precision                              , intent(in   ), optional :: massThreshold
    integer                                                                 :: i               , indexMinimum, &
         &                                                                     indexMaximum
    double precision                                                        :: factor

    call fileTableIndexRange(table,indexMetallicity,massInitial,massThreshold,indexMinimum,indexMaximum)
    if (massEvaluate >= table%mass(indexMaximum,indexMetallicity)) then
       ! Extend beyond the most massive star with a fixed value per unit initial mass.
       value =+table%value       (indexMaximum,indexMetallicity) &
            & *      massEvaluate                                &
            & /table%mass        (indexMaximum,indexMetallicity)
       return
    end if
    do i=indexMinimum,indexMaximum-2
       if (massEvaluate < table%mass(i+1,indexMetallicity)) exit
    end do
    factor=+(      massEvaluate                      -table%mass(i,indexMetallicity)) &
         & /(table%mass        (i+1,indexMetallicity)-table%mass(i,indexMetallicity))
    value =+(1.0d0-factor)*table%value(i  ,indexMetallicity) &
         & +       factor *table%value(i+1,indexMetallicity)
    return
  end function fileTableInterpolateMass

  subroutine fileValidateMetallicities(fileName,label,metallicities)
    !!{RST
    Report a property which is tabulated at only a single metallicity.

    Such a property cannot be interpolated: the tabulated points are collinear in the (mass, metallicity) plane,
    so interpolation there is not defined and any request at another metallicity is an extrapolation. The
    failure is silent---the interpolation returns a value, just not a meaningful one, and has been observed to
    return an order of magnitude too much, or zero, or values which send the integration over the initial mass
    function into an unbounded subdivision---so it is caught here rather than left to be discovered downstream.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    implicit none
    type            (varying_string), intent(in   )               :: fileName
    character       (len=*         ), intent(in   )               :: label
    double precision                , intent(in   ), dimension(:) :: metallicities
    character       (len=24        )                              :: metallicityLabel

    if (size(metallicities) <= 1                       ) return
    if (maxval(metallicities) > minval(metallicities)  ) return
    write (metallicityLabel,'(e12.6)') metallicities(1)
    call Error_Report(                                                                                                  &
         &            'stellar properties file "'//char(fileName)//'" tabulates '//trim(label)//                        &
         &            ' at a single metallicity (Z='//trim(adjustl(metallicityLabel))//')'//char(10)//                  &
         &            '   => such a property can not be interpolated in the (mass, metallicity) plane, and a request'// &
         &            ' at any other metallicity would be extrapolated, returning a value which is not meaningful'   // &
         &            char(10)//                                                                                        &
         &            '   => if the property genuinely does not depend on metallicity, tabulate it at two'           // &
         &            ' metallicities which bracket the range of interest, giving it the same value at both; the'    // &
         &            ' Sukhbold et al. (2016) files are written this way'                                           // &
         &            {introspection:location}                                                                          &
         &           )
    return
  end subroutine fileValidateMetallicities

  double precision function fileMassInitial(self,lifetime,metallicity)
    !!{RST
    Return the initial mass of a star of given ``lifetime`` and ``metallicity``.
    !!}
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
         &                                   dataStatic         =.true.                                , &
         &                                   numberComputePoints=     3                                  &
         &                                  )
    return
  end function fileMassInitial

  double precision function fileLifetime(self,massInitial,metallicity)
    !!{RST
    Return the lifetime of a star (in Gyr) given an ``massInitial`` and ``metallicity``.
    !!}
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
         &                                reset=self%interpolationResetLifetime    , &
         &                                dataStatic=.true.                          &
         &                               )
    return
  end function fileLifetime

  double precision function fileMassEjected(self,massInitial,metallicity)
    !!{RST
    Return the mass ejected during the lifetime of a star of given ``massInitial`` and ``metallicity``.
    !!}
    implicit none
    class           (stellarAstrophysicsFile), intent(inout) :: self
    double precision                         , intent(in   ) :: massInitial, metallicity

    call self%read()
    fileMassEjected=max(fileTableInterpolate(self%tableMassEjected,massInitial,metallicity,self%initialMassForSupernovaeTypeII),0.0d0)
    return
  end function fileMassEjected

  double precision function fileMassYield(self,massInitial,metallicity,atomIndex)
    !!{RST
    Return the mass of metals (or, if ``atomIndex`` is given, of that element) yielded by a star of given ``massInitial`` and
    ``metallicity``.
    !!}
    implicit none
    class           (stellarAstrophysicsFile), intent(inout)           :: self
    double precision                         , intent(in   )           :: massInitial , metallicity
    integer                                  , intent(in   ), optional :: atomIndex
    integer                                                            :: elementIndex

    call self%read()
    if (present(atomIndex)) then
       ! Compute the element mass yield. Elements for which this compilation provides no yields at all are mapped to an index
       ! of -1, and have zero yield.
       elementIndex=self%atomIndexMap(atomIndex)
       if (elementIndex > 0) then
          fileMassYield=max(fileTableInterpolate(self%tableYieldElement(elementIndex),massInitial,metallicity,self%initialMassForSupernovaeTypeII),0.0d0)
       else
          fileMassYield=0.0d0
       end if
    else
       ! Compute the metal mass yield.
       fileMassYield   =max(fileTableInterpolate(self%tableYieldMetals               ,massInitial,metallicity,self%initialMassForSupernovaeTypeII),0.0d0)
    end if
    return
  end function fileMassYield

  double precision function fileMassInitialSupernovaeTypeII(self)
    !!{RST
    Return the minimum initial mass of a star which ends its life as a Type II supernova.
    !!}
    implicit none
    class(stellarAstrophysicsFile), intent(inout) :: self

    fileMassInitialSupernovaeTypeII=self%initialMassForSupernovaeTypeII
    return
  end function fileMassInitialSupernovaeTypeII
