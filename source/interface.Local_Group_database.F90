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

!!{
Contains a module which interfaces with the Local Group database.
!!}

module Interface_Local_Group_DB
  !!{
  Interfaces with the Local Group database.
  !!}
  use :: FoX_DOM         , only : node
  use :: IO_XML          , only : xmlNodeList
  use :: Resource_Manager, only : resourceManager
  private
  public  :: localGroupDB, vector3D

  !![
  <generic identifier="Type">
   <instance label="Double" intrinsic="double precision"                />
   <instance label="VarStr" intrinsic="type            (varying_string)"/>
   <instance label="Vector" intrinsic="type            (vector3D      )"/>
  </generic>
  !!]

  !![
  <enumeration>
   <name>comparison</name>
   <description>Comparison operators.</description>
   <visibility>public</visibility>
   <entry label="equals"     />
   <entry label="greaterThan"/>
   <entry label="lessThan"   />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>setOperator</name>
   <description>Set operators.</description>
   <visibility>public</visibility>
   <entry label="intersection"/>
   <entry label="union"/>
   <entry label="relativeComplement"/>
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>attribute</name>
   <description>Attribute types.</description>
   <visibility>public</visibility>
   <entry label="value"      />
   <entry label="uncertainty"/>
  </enumeration>
  !!]

  type :: vector3D
     !!{
     Vector type.
     !!}
     double precision, dimension(3) :: x
   contains
     !![
     <methods>
      <method method="operator(==)"   description="Test equality of two 3D vectors."         />
      <method method="operator(&lt;)" description="Less than operator for two 3D vectors."   />
      <method method="operator(&gt;)" description="Greater than operator for two 3D vectors."/>
     </methods>
     !!]
     procedure ::                 vector3DEquals
     procedure ::                 vector3DComparisonUnimplemented
     generic   :: operator(==) => vector3DEquals
     generic   :: operator(<)  => vector3DComparisonUnimplemented
     generic   :: operator(>)  => vector3DComparisonUnimplemented
  end type vector3D

  type :: documentWrapper
     !!{
     Wrapper class for managing XML documents.
     !!}
     private
     type(node), pointer :: document => null()
   contains
     final :: documentWrapperDestructor
  end type documentWrapper
  
  type :: localGroupDB
     !!{
     Local Group database class.
     !!}
     private
     type   (documentWrapper), pointer                   :: database        => null()
     type   (xmlNodeList    ), allocatable, dimension(:) :: galaxies
     logical                 , allocatable, dimension(:) :: selected
     type   (resourceManager)                            :: databaseManager
   contains
     !![
     <methods>
       <method method="getProperty" description="Return an array of values of the named property for the current selection."                />
       <method method="select"      description="Select all galaxies in the current selection where the named property has the given value."/>
       <method method="selectAll"   description="Select all galaxies in the database."                                                      />
       <method method="update"      description="Update the database."                                                                      />
     </methods>
     !!]
     procedure :: getProperty{Type¦label} => localGroupDBGetProperty{Type¦label}
     generic   :: getProperty             => getProperty{Type¦label}
     procedure :: select{Type¦label}      => localGroupDBSelect{Type¦label}
     generic   :: select                  => select{Type¦label}
     procedure :: selectAll               => localGroupDBSelectAll
     procedure :: update                  => localGroupDBUpdate
  end type localGroupDB

  interface localGroupDB
     !!{
     Constructors for the Local Group database class.
     !!}
     module procedure localGroupDBConstructorInternal
  end interface localGroupDB

contains

  function localGroupDBConstructorInternal() result(self)
    !!{
    Constructor for the Local Group database class.
    !!}
    use :: Input_Paths       , only : inputPath                   , pathTypeDataStatic
    use :: IO_XML            , only : XML_Get_Elements_By_Tag_Name, XML_Get_First_Element_By_Tag_Name, XML_Parse
    use :: ISO_Varying_String, only : char                        , varying_string                   , operator(//)
    implicit none
    type (localGroupDB  )          :: self
    type (varying_string)          :: fileName
    class(*             ), pointer :: dummyPointer_

    fileName=inputPath(pathTypeDataStatic)//"observations/localGroup/localGroupSatellites.xml"
    allocate(self%database)
    !$omp critical (FoX_DOM_Access)
    self%database%document => XML_Parse(char(fileName))
    call XML_Get_Elements_By_Tag_Name(XML_Get_First_Element_By_Tag_Name(self%database%document,'galaxies'),'galaxy',self%galaxies)
    !$omp end critical (FoX_DOM_Access)
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_ => self%database
    self%databaseManager=resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    allocate(self%selected(0:size(self%galaxies)-1))
    self%selected=.false.
    return
  end function localGroupDBConstructorInternal

  subroutine documentWrapperDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily documentWrapper} object.
    !!}
    use :: FoX_DOM, only : destroy
    implicit none
    type(documentWrapper), intent(inout) :: self
    
    if (associated(self%document)) then
       !$omp critical (FoX_DOM_Access)
       call destroy(self%document)
       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end subroutine documentWrapperDestructor

  subroutine localGroupDBGetProperty{Type¦label}(self,name,property,isPresent,attribute)
    !!{
    Get a named text property from the Local Group database.
    !!}
    use                      :: FoX_DOM           , only : getAttributeNode            , getTextContent , hasAttribute, setAttribute, &
         &                                                 extractDataContent
    use                      :: Error             , only : Error_Report
    use                      :: IO_XML            , only : XML_Get_Elements_By_Tag_Name
    {Type¦match¦^VarStr$¦use :: ISO_Varying_String, only : varying_string              , assignment(=)¦}
    implicit none
    class           (localGroupDB            ), intent(inout)                                      :: self
    character       (len=*                   ), intent(in   )                                      :: name
    {Type¦intrinsic}                          , intent(inout), allocatable, dimension(:)           :: property
    logical                                   , intent(inout), allocatable, dimension(:), optional :: isPresent
    type            (enumerationAttributeType), intent(in   )                           , optional :: attribute
    type            (node                    )               , pointer                             :: galaxy          , propertyElement, &
         &                                                                                            attributeNode
    type            (xmlNodeList             )               , allocatable, dimension(:)           :: properties
    integer                                                                                        :: i               , countGalaxies
    logical                                                                                        :: propertyFound
    double precision                                                                               :: uncertaintyLower, uncertaintyUpper, &
         &                                                                                            uncertainty
    character       (len=64                  )                                                     :: textContent
    !![
    <optionalArgument name="attribute" defaultsTo="attributeValue" />
    !!]
    
    ! Iterate over galaxies to count how many have this property.
    countGalaxies=0
    do i=0,size(self%galaxies)-1
       if (.not.self%selected(i)) cycle
       galaxy => self%galaxies(i)%element
       if (present(isPresent)) then
          countGalaxies=countGalaxies+1
       else if (attribute_ == attributeValue .and. hasAttribute(galaxy,name)) then
          countGalaxies=countGalaxies+1
       else
          call XML_Get_Elements_By_Tag_Name(galaxy,name,properties)
          if      (size(properties) > 1) then
             call Error_Report('galaxy has multiple entries for named property'//{introspection:location})
          else if (size(properties) == 1) then
             propertyElement => properties(0)%element
             if     (                                                                                                                                                 &
                  &   (attribute_ == attributeValue       .and. hasAttribute(propertyElement,'value'         )                                                      ) &
                  &  .or.                                                                                                                                             &
                  &   (attribute_ == attributeUncertainty .and. hasAttribute(propertyElement,'uncertainty'   )                                                      ) &
                  &  .or.                                                                                                                                             &
                  &   (attribute_ == attributeUncertainty .and. hasAttribute(propertyElement,'uncertaintyLow') .and. hasAttribute(propertyElement,'uncertaintyHigh')) &
                  & ) then
                countGalaxies=countGalaxies+1
                if (attribute_ == attributeUncertainty .and. .not.hasAttribute(propertyElement,'uncertainty')) then
                   ! Compute an uncertainty from the lower and upper bounds.
                   propertyElement => properties(0)%element
                   attributeNode   => getAttributeNode(propertyElement,'uncertaintyLower')
                   call extractDataContent(attributeNode,uncertaintyLower)
                   attributeNode   => getAttributeNode(propertyElement,'uncertaintyUpper')
                   call extractDataContent(attributeNode,uncertaintyUpper)
                   uncertainty=0.5d0*(uncertaintyLower+uncertaintyUpper)
                   write (textContent,'(f16.12)') uncertainty
                   call setAttribute(propertyElement,"uncertainty",trim(adjustl(textContent)))
                end if
             end if
          end if
       end if
    end do
    allocate(property(countGalaxies))
    if (present(isPresent)) then
       allocate(isPresent(countGalaxies))
       isPresent=.true.
    end if
    ! Populate the array.
    countGalaxies=0
    do i=0,size(self%galaxies)-1
       if (.not.self%selected(i)) cycle
       countGalaxies =  countGalaxies+1
       propertyFound =  .false.
       galaxy        => self%galaxies(i)%element
       if (attribute_ == attributeValue .and. hasAttribute(galaxy,name)) then
          propertyFound =  .true.
          attributeNode => getAttributeNode(galaxy,name)
       else
          call XML_Get_Elements_By_Tag_Name(galaxy,name,properties)
          if (size(properties) == 1) then
             propertyElement => properties(0)%element
             if     (                                                                                           &
                  &   (attribute_ == attributeValue       .and. hasAttribute(propertyElement,'value'         )) &
                  &  .or.                                                                                       &
                  &   (attribute_ == attributeUncertainty .and. hasAttribute(propertyElement,'uncertainty'   )) &
                  & ) then
                propertyFound=.true.
                if       (attribute_ == attributeValue     ) then
                   attributeNode => getAttributeNode(propertyElement,'value'      )
                else  if (attribute_ == attributeUncertainty) then
                   attributeNode => getAttributeNode(propertyElement,'uncertainty')
                end if
             end if
          end if
       end if
       if (propertyFound) then
          {Type¦match¦^Double$¦call extractDataContent(attributeNode,property(countGalaxies))¦}
          {Type¦match¦^Vector¦call extractDataContent(attributeNode,property(countGalaxies)%x)¦}
          {Type¦match¦^VarStr$¦property(countGalaxies)=getTextContent(attributeNode)¦}
       else
          if (present(isPresent)) then
             isPresent(countGalaxies)=.false.
             {Type¦match¦^Double$¦property(countGalaxies)=0.0d0¦}
             {Type¦match¦^Vector¦property(countGalaxies)%x=0.0d0¦}
             {Type¦match¦^VarStr$¦property(countGalaxies)='not present'¦}
          else
             attributeNode => getAttributeNode(galaxy,'name')
             call Error_Report('property "'//name//'" is not present in galaxy "'//getTextContent(attributeNode)//'"'//{introspection:location})
          end if
       end if       
    end do
    return
  end subroutine localGroupDBGetProperty{Type¦label}

  subroutine localGroupDBSelectAll(self)
    !!{
    Select all galaxies in the database.
    !!}
    implicit none
    class(localGroupDB), intent(inout) :: self

    self%selected=.true.
    return
  end subroutine localGroupDBSelectAll

  subroutine localGroupDBSelect{Type¦label}(self,name,value,comparison,setOperator)
    !!{
    Impose a selection on the database.
    !!}
    use                      :: FoX_DOM           , only : getAttributeNode, getTextContent
    use                      :: Error             , only : Error_Report
    use                      :: ISO_Varying_String, only : varying_string
    {Type¦match¦^VarStr$¦use :: ISO_Varying_String, only : operator(<)     , operator(>)   , operator(==)¦}
    implicit none
    class           (localGroupDB              ), intent(inout)               :: self
    character       (len=*                     ), intent(in   )               :: name
    {Type¦intrinsic}                            , intent(in)                  :: value
    type            (enumerationComparisonType ), intent(in   )               :: comparison
    type            (enumerationSetOperatorType), intent(in   )               :: setOperator
    {Type¦intrinsic}                            , allocatable  , dimension(:) :: values
    logical                                     , allocatable  , dimension(:) :: selectedCurrent , isPresent
    type            (node                      )               , pointer      :: galaxy          , attribute
    integer                                                                   :: i
    logical                                                                   :: comparisonResult

    allocate(selectedCurrent(0:size(self%galaxies)-1))
    selectedCurrent=self%selected
    self%selected=.true.
    call self%getProperty(name,values,isPresent)
    do i=0,size(self%galaxies)-1
       if (.not.selectedCurrent(i) .and. (setOperator == setOperatorIntersection .or.  setOperator == setOperatorRelativeComplement)) cycle
       if (.not.isPresent(i+1)) then
          galaxy    => self%galaxies(i)%element
          attribute => getAttributeNode(galaxy,'name')
         call Error_Report('property "'//name//'" is not present in selected galaxy "'//getTextContent(attribute)//'"'//{introspection:location})
       end if
       select case (comparison%ID)
       case (comparisonEquals     %ID)
          comparisonResult=values(i+1) == value
       case (comparisonLessThan   %ID)
          comparisonResult=values(i+1) <  value
       case (comparisonGreaterThan%ID)
          comparisonResult=values(i+1) >  value
       case default
          comparisonResult=.false.
          call Error_Report('unknown comparison operator'//{introspection:location})
       end select
       select case (setOperator%ID)
       case (setOperatorIntersection      %ID)
          selectedCurrent(i)=selectedCurrent(i) .and.      comparisonResult
       case (setOperatorUnion             %ID)
          selectedCurrent(i)=selectedCurrent(i) .or.       comparisonResult
       case (setOperatorRelativeComplement%ID)
          selectedCurrent(i)=selectedCurrent(i) .and. .not.comparisonResult
       case default
          call Error_Report('unknown set operator'       //{introspection:location})
       end select
    end do
    self%selected=selectedCurrent
    return
  end subroutine localGroupDBSelect{Type¦label}

  subroutine localGroupDBUpdate(self)
    !!{
    Update the database.
    !!}
    use :: FoX_DOM                         , only : appendChild                 , createElementNS    , setAttribute    , getAttributeNode, &
         &                                          getNamespaceURI             , getTextContent     , hasAttribute    , serialize       , &
         &                                          extractDataContent
    use :: Error                           , only : Error_Report
    use :: Input_Paths                     , only : inputPath                   , pathTypeDataStatic
    use :: IO_XML                          , only : XML_Get_Elements_By_Tag_Name
    use :: ISO_Varying_String              , only : char
    use :: Numerical_Constants_Astronomical, only : arcminutesToDegrees         , arcsecondsToDegrees, degreesToRadians, hoursToDegrees  , &
          &                                         minutesToDegrees            , secondsToDegrees
    implicit none
    class           (localGroupDB  ), intent(inout)             :: self
    type            (node          ), pointer                   :: galaxy                        , attribute                       , newNode
    type            (xmlNodeList   ), dimension(:), allocatable :: propertyList1                 , propertyList2                   , propertyList3                   , propertyList4
    character       (len=32        ), dimension(4)              :: uncertainties
    double precision                , dimension(3)              :: positionMilkyWay              , positionM31                     , position                        , uncertaintyPosition  , &
         &                                                         uncertaintyPositionMilkyWay   , uncertaintyPositionM31
    integer                                                     :: i                             , j                               , indexMilkyWay                   , indexM31
    logical                                                     :: hasSexagesimal                , hasDecimal                      , hasHeliocentric                 , hasModulus           , &
         &                                                         hasDistance                   , hasDeclination                  , hasRightAscension               , hasPosition
    double precision                                            :: declinationSexagesimalDegrees , declinationSexagesimalArcMinutes, declinationSexagesimalArcSeconds, declinationDecimal   , &
         &                                                         rightAscensionSexagesimalHours, rightAscensionSexagesimalMinutes, rightAscensionSexagesimalSeconds, rightAscensionDecimal, &
         &                                                         distanceHeliocentric          , distanceModulus                 , uncertainty1                    , uncertainty2         , &
         &                                                         positionHeliocentricX         , positionHeliocentricY           , positionHeliocentricZ           , distance
    character       (len=64        )                            :: textContent

    ! Set names of uncertainties.
    uncertainties(1)='uncertainty'
    uncertainties(2)='uncertaintyLow'
    uncertainties(3)='uncertaintyHigh'
    uncertainties(4)='uncertaintySystematic'
    ! Identify Milky Way and M31.
    indexMilkyWay=-huge(0)
    indexM31     =-huge(0)
    do i=0,size(self%galaxies)-1
       galaxy    => self%galaxies(i)%element
       attribute => getAttributeNode(galaxy,'name')
       if (getTextContent(attribute) == "The Galaxy") indexMilkyWay=i
       if (getTextContent(attribute) == "Andromeda" ) indexM31     =i
    end do
    if (indexMilkyWay < 0) call Error_Report('unable to find Milky Way in the database'//{introspection:location})
    if (indexM31      < 0) call Error_Report('unable to find M31 in the database'      //{introspection:location})
    ! Iterate over galaxies.
    do i=0,size(self%galaxies)-1
       galaxy => self%galaxies(i)%element
       ! Look for declinations.
       call XML_Get_Elements_By_Tag_Name(galaxy,'declinationSexagesimal',propertyList1)
       hasSexagesimal=size(propertyList1) == 1
       call XML_Get_Elements_By_Tag_Name(galaxy,'declinationDecimal'    ,propertyList2)
       hasDecimal    =size(propertyList2) == 1
       if (hasSexagesimal .and. .not. hasDecimal) then
          ! Convert sexagesimal declination to decimal.
          attribute => getAttributeNode(propertyList1(0)%element,'degrees'  )
          call extractDataContent(attribute,declinationSexagesimalDegrees   )
          attribute => getAttributeNode(propertyList1(0)%element,'arcminutes')
          call extractDataContent(attribute,declinationSexagesimalArcMinutes)
          attribute => getAttributeNode(propertyList1(0)%element,'arcseconds')
          call extractDataContent(attribute,declinationSexagesimalArcSeconds)
          declinationDecimal=+     declinationSexagesimalDegrees                                                       &
               &             +sign(declinationSexagesimalArcMinutes,declinationSexagesimalDegrees)*arcminutesToDegrees &
               &             +sign(declinationSexagesimalArcSeconds,declinationSexagesimalDegrees)*arcsecondsToDegrees
          newNode => createElementNS(self%database%document,getNamespaceURI(self%database%document),'declinationDecimal')
          write (textContent,'(f16.12)') declinationDecimal
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          newNode => appendChild(galaxy,newNode)
       else if (.not.hasSexagesimal .and. hasDecimal) then
          ! Convert decimal declination to sexagesimal.
          attribute => getAttributeNode(propertyList2(0)%element,'value')
          call extractDataContent(attribute,declinationDecimal)
          declinationSexagesimalDegrees   =int( declinationDecimal                                                   )
          declinationSexagesimalArcminutes=int((declinationDecimal-declinationSexagesimalDegrees)/arcminutesToDegrees)
          declinationSexagesimalArcseconds=    (declinationDecimal-declinationSexagesimalDegrees-declinationSexagesimalArcminutes*arcminutesToDegrees)/arcsecondsToDegrees
          declinationSexagesimalArcminutes=abs(declinationSexagesimalArcminutes)
          declinationSexagesimalArcseconds=abs(declinationSexagesimalArcseconds)
          newNode => createElementNS(self%database%document,getNamespaceURI(self%database%document),'declinationSexagesimal')
          write (textContent,'(f16.12)') declinationSexagesimalDegrees
          call setAttribute(newNode,"degrees"   ,trim(adjustl(textContent)))
          write (textContent,'(f16.12)') declinationSexagesimalArcMinutes
          call setAttribute(newNode,"arcminutes",trim(adjustl(textContent)))
          write (textContent,'(f16.12)') declinationSexagesimalArcSeconds
          call setAttribute(newNode,"arcseconds",trim(adjustl(textContent)))
          newNode => appendChild(galaxy,newNode)
       end if
       ! Look for right ascensions.
       call XML_Get_Elements_By_Tag_Name(galaxy,'rightAscensionSexagesimal',propertyList1)
       hasSexagesimal=size(propertyList1) == 1
       call XML_Get_Elements_By_Tag_Name(galaxy,'rightAscensionDecimal'    ,propertyList2)
       hasDecimal=size(propertyList2) == 1
       if (hasSexagesimal .and. .not. hasDecimal) then
          ! Convert sexagesimal right ascension to decimal.
          attribute => getAttributeNode(propertyList1(0)%element,'hours'  )
          call extractDataContent(attribute,rightAscensionSexagesimalHours  )
          attribute => getAttributeNode(propertyList1(0)%element,'minutes')
          call extractDataContent(attribute,rightAscensionSexagesimalMinutes)
          attribute => getAttributeNode(propertyList1(0)%element,'seconds')
          call extractDataContent(attribute,rightAscensionSexagesimalSeconds)
          rightAscensionDecimal=+rightAscensionSexagesimalHours  *hoursToDegrees   &
               &                +rightAscensionSexagesimalMinutes*minutesToDegrees &
               &                +rightAscensionSexagesimalSeconds*secondsToDegrees
          newNode => createElementNS(self%database%document,getNamespaceURI(self%database%document),'rightAscensionDecimal')
          write (textContent,'(f16.12)') rightAscensionDecimal
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          newNode => appendChild(galaxy,newNode)
       else if (.not.hasSexagesimal .and. hasDecimal) then
          ! Convert decimal right ascension to sexagesimal.
          attribute => getAttributeNode(propertyList2(0)%element,'value')
          call extractDataContent(attribute,rightAscensionDecimal)
          rightAscensionSexagesimalHours  =int( rightAscensionDecimal/hoursToDegrees                                                 )
          rightAscensionSexagesimalMinutes=int((rightAscensionDecimal-rightAscensionSexagesimalHours*hoursToDegrees)/minutesToDegrees)
          rightAscensionSexagesimalSeconds=(rightAscensionDecimal-rightAscensionSexagesimalHours*hoursToDegrees-rightAscensionSexagesimalMinutes*minutesToDegrees)/secondsToDegrees
          newNode => createElementNS(self%database%document,getNamespaceURI(self%database%document),'rightAscensionSexagesimal')
          write (textContent,'(f16.12)') rightAscensionSexagesimalHours
          call setAttribute(newNode,"hours"  ,trim(adjustl(textContent)))
          write (textContent,'(f16.12)') rightAscensionSexagesimalMinutes
          call setAttribute(newNode,"minutes",trim(adjustl(textContent)))
          write (textContent,'(f16.12)') rightAscensionSexagesimalSeconds
          call setAttribute(newNode,"seconds",trim(adjustl(textContent)))
          newNode => appendChild(galaxy,newNode)
       end if
       ! Look for distances
       call XML_Get_Elements_By_Tag_Name(galaxy,'distanceHeliocentric',propertyList1)
       hasHeliocentric=size(propertyList1) == 1
       call XML_Get_Elements_By_Tag_Name(galaxy,'distanceModulus'     ,propertyList2)
       hasModulus=size(propertyList2) == 1
       if (hasHeliocentric .and. .not. hasModulus) then
          ! Convert heliocentric distance to a distance modulus.
          attribute => getAttributeNode(propertyList1(0)%element,'value')
          call extractDataContent(attribute,distanceHeliocentric)
          distanceModulus=25.0d0+5.0d0*log10(distanceHeliocentric)
          newNode => createElementNS(self%database%document,getNamespaceURI(self%database%document),'distanceModulus')
          write (textContent,'(f16.12)') distanceModulus
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          do j=1,size(uncertainties)
             if (hasAttribute(propertyList1(0)%element,trim(uncertainties(j)))) then
                attribute => getAttributeNode(propertyList1(0)%element,trim(uncertainties(j)))
                call extractDataContent(attribute,uncertainty1)
                uncertainty2=5.0d0/log(10.0d0)*uncertainty1/distanceHeliocentric
                write (textContent,'(f16.12)') uncertainty2
                call setAttribute(newNode,trim(uncertainties(j)),trim(adjustl(textContent)))
             end if
          end do
          newNode => appendChild(galaxy,newNode)
       else if (.not.hasHeliocentric .and. hasModulus) then
          ! Convert distance modulus to a heliocentric distance.
          attribute => getAttributeNode(propertyList2(0)%element,'value')
          call extractDataContent(attribute,distanceModulus)
          distanceHeliocentric=10.0d0**((distanceModulus-25.0d0)/5.0d0)
          newNode => createElementNS(self%database%document,getNamespaceURI(self%database%document),'distanceHeliocentric')
          write (textContent,'(f16.12)') distanceHeliocentric
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          do j=1,size(uncertainties)
             if (hasAttribute(propertyList2(0)%element,trim(uncertainties(j)))) then
                attribute => getAttributeNode(propertyList2(0)%element,trim(uncertainties(j)))
                call extractDataContent(attribute,uncertainty1)
                uncertainty2=distanceHeliocentric*log(10.0d0)/5.0d0
                write (textContent,'(f16.12)') uncertainty2
                call setAttribute(newNode,trim(uncertainties(j)),trim(adjustl(textContent)))
             end if
          end do
          newNode => appendChild(galaxy,newNode)
       end if
       ! Compute Cartesian heliocentric coordinates.
       call XML_Get_Elements_By_Tag_Name(galaxy,'distanceHeliocentric'         ,propertyList1)
       hasDistance=size(propertyList1) == 1
       call XML_Get_Elements_By_Tag_Name(galaxy,'rightAscensionDecimal'        ,propertyList2)
       hasRightAscension=size(propertyList2) == 1
       call XML_Get_Elements_By_Tag_Name(galaxy,'declinationDecimal'           ,propertyList3)
       hasDeclination=size(propertyList3) == 1
       call XML_Get_Elements_By_Tag_Name(galaxy,'positionHeliocentricCartesian',propertyList4)
       hasPosition=size(propertyList4) == 1
       if (hasDistance .and. hasRightAscension .and. hasDeclination .and. .not. hasPosition) then
          attribute => getAttributeNode(propertyList1(0)%element,'value')
          call extractDataContent(attribute,distanceHeliocentric )
          attribute => getAttributeNode(propertyList2(0)%element,'value')
          call extractDataContent(attribute,rightAscensionDecimal)
          attribute => getAttributeNode(propertyList3(0)%element,'value')
          call extractDataContent(attribute,declinationDecimal   )
          positionHeliocentricX=distanceHeliocentric*cos(declinationDecimal*degreesToRadians)*cos(rightAscensionDecimal*degreesToRadians)
          positionHeliocentricY=distanceHeliocentric*cos(declinationDecimal*degreesToRadians)*sin(rightAscensionDecimal*degreesToRadians)
          positionHeliocentricZ=distanceHeliocentric*sin(declinationDecimal*degreesToRadians)
          newNode => createElementNS(self%database%document,getNamespaceURI(self%database%document),'positionHeliocentricCartesian')
          write (textContent,'(e16.8,1x,e16.8,1x,e16.8)') positionHeliocentricX,positionHeliocentricY,positionHeliocentricZ
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          do j=1,size(uncertainties)
             if (hasAttribute(propertyList1(0)%element,trim(uncertainties(j)))) then
                attribute => getAttributeNode(propertyList1(0)%element,trim(uncertainties(j)))
                call extractDataContent(attribute,uncertainty1)
                positionHeliocentricX=abs(uncertainty1*cos(declinationDecimal*degreesToRadians)*cos(rightAscensionDecimal*degreesToRadians))
                positionHeliocentricY=abs(uncertainty1*cos(declinationDecimal*degreesToRadians)*sin(rightAscensionDecimal*degreesToRadians))
                positionHeliocentricZ=abs(uncertainty1*sin(declinationDecimal*degreesToRadians))
                write (textContent,'(e16.8,1x,e16.8,1x,e16.8)') positionHeliocentricX,positionHeliocentricY,positionHeliocentricZ
                call setAttribute(newNode,trim(uncertainties(j)),trim(adjustl(textContent)))
             end if
          end do
          newNode => appendChild(galaxy,newNode)
       end if
    end do
    ! Get heliocentric Cartesian positions of Milky Way and M31.
    galaxy   => self%galaxies(indexMilkyWay)%element
    call XML_Get_Elements_By_Tag_Name(galaxy,'positionHeliocentricCartesian',propertyList1)
    attribute => getAttributeNode(propertyList1(0)%element,'value'      )
    call extractDataContent(attribute,positionMilkyWay           )
    attribute => getAttributeNode(propertyList1(0)%element,'uncertainty')
    call extractDataContent(attribute,uncertaintyPositionMilkyWay)
    galaxy  => self%galaxies(indexM31     )%element
    call XML_Get_Elements_By_Tag_Name(galaxy,'positionHeliocentricCartesian',propertyList1)
    attribute => getAttributeNode(propertyList1(0)%element,'value'      )
    call extractDataContent(attribute,positionM31)
    attribute => getAttributeNode(propertyList1(0)%element,'uncertainty')
    call extractDataContent(attribute,uncertaintyPositionM31)
    ! Compute galactocentric distances.
    do i=0,size(self%galaxies)-1
       ! Milky Way.
       galaxy => self%galaxies(i)%element
       call XML_Get_Elements_By_Tag_Name(galaxy,'positionHeliocentricCartesian',propertyList1)
       hasPosition=size(propertyList1) == 1
       call XML_Get_Elements_By_Tag_Name(galaxy,'distanceMilkyWay'             ,propertyList2)
       hasDistance=size(propertyList2) == 1
       if (hasPosition .and. .not.hasDistance) then
          attribute => getAttributeNode(propertyList1(0)%element,'value')
          call extractDataContent(attribute,position)
          distance=sqrt(sum((position-positionMilkyWay)**2))
          newNode => createElementNS(self%database%document,getNamespaceURI(self%database%document),'distanceMilkyWay')
          write (textContent,'(e16.8)') distance
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          do j=1,size(uncertainties)
             if (hasAttribute(propertyList1(0)%element,trim(uncertainties(j)))) then
                attribute => getAttributeNode(propertyList1(0)%element,trim(uncertainties(j)))
                call extractDataContent(attribute,uncertaintyPosition)
                if (distance > 0.0d0) then
                   uncertainty2=sqrt(sum((position-positionMilkyWay)**2*(uncertaintyPosition**2+uncertaintyPositionMilkyWay**2)))/distance
                else
                   uncertainty2=0.0d0
                end if
                write (textContent,'(e16.8)') uncertainty2
                call setAttribute(newNode,trim(uncertainties(j)),trim(adjustl(textContent)))
             end if
          end do
          newNode => appendChild(galaxy,newNode)
       end if
       ! M31.
       galaxy => self%galaxies(i)%element
       call XML_Get_Elements_By_Tag_Name(galaxy,'positionHeliocentricCartesian',propertyList1)
       hasPosition=size(propertyList1) == 1
       call XML_Get_Elements_By_Tag_Name(galaxy,'distanceM31'                  ,propertyList2)
       hasDistance=size(propertyList2) == 1
       if (hasPosition .and. .not.hasDistance) then
          attribute => getAttributeNode(propertyList1(0)%element,'value')
          call extractDataContent(attribute,position)
          distance=sqrt(sum((position-positionM31)**2))
          newNode => createElementNS(self%database%document,getNamespaceURI(self%database%document),'distanceM31')
          write (textContent,'(e16.8)') distance
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          do j=1,size(uncertainties)
             if (hasAttribute(propertyList1(0)%element,trim(uncertainties(j)))) then
                attribute => getAttributeNode(propertyList1(0)%element,trim(uncertainties(j)))
                attribute => getAttributeNode(propertyList1(0)%element,trim(uncertainties(j)))
                call extractDataContent(attribute,uncertaintyPosition)
                if (distance > 0.0d0) then
                   uncertainty2=sqrt(sum((position-positionM31)**2*(uncertaintyPosition**2+uncertaintyPositionM31**2)))/distance
                else
                   uncertainty2=0.0d0
                end if
                write (textContent,'(e16.8)') uncertainty2
                call setAttribute(newNode,trim(uncertainties(j)),trim(adjustl(textContent)))
             end if
          end do
          newNode => appendChild(galaxy,newNode)
       end if
    end do
    ! Write out the updated database.
    call serialize(self%database%document,char(inputPath(pathTypeDataStatic))//"observations/localGroup/localGroupSatellites.xml")
    return
  end subroutine localGroupDBUpdate

  logical function vector3DEquals(self,other)
    !!{
    Equality comparison operator for two 3D vectors.
    !!}
    implicit none
    class(vector3D), intent(in   ) :: self, other

    vector3DEquals=all(self%x == other%x)
    return
  end function vector3DEquals

  logical function vector3DComparisonUnimplemented(self,other)
    !!{
    Unimplemented comparison operators for 3D vectors.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(vector3D), intent(in   ) :: self, other
    !$GLC attributes unused :: self, other

    vector3DComparisonUnimplemented=.false.
    call Error_Report('comparison operator is unimplemented'//{introspection:location})
    return
  end function vector3DComparisonUnimplemented

end module Interface_Local_Group_DB
