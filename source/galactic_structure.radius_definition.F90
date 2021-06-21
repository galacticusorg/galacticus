!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which provides parsing of radii definitions used in output specifiers.
!!}

module Galactic_Structure_Radii_Definitions
  !!{
  Provides parsing of radii definitions used in output specifiers.
  !!}
  use :: ISO_Varying_String, only : varying_string
  private
  public :: radiusSpecifier, Galactic_Structure_Radii_Definition_Decode

  type radiusSpecifier
     !!{
     Type used for specifying radii definitions used in output specifiers.
     !!}
     type            (varying_string) :: name
     integer                          :: component            , direction    , integralWeightBy, &
          &                              integralWeightByIndex, mass         , type            , &
          &                              weightBy             , weightByIndex
     double precision                 :: fraction             , value
  end type radiusSpecifier

  !![
  <enumeration>
   <name>radiusType</name>
   <description>Used to specify radii types used in output specifiers.</description>
   <encodeFunction>yes</encodeFunction>
   <visibility>public</visibility>
   <validator>yes</validator>
   <entry label="radius"                />
   <entry label="virialRadius"          />
   <entry label="darkMatterScaleRadius" />
   <entry label="diskRadius"            />
   <entry label="spheroidRadius"        />
   <entry label="diskHalfMassRadius"    />
   <entry label="spheroidHalfMassRadius"/>
   <entry label="galacticMassFraction"  />
   <entry label="galacticLightFraction" />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>direction</name>
   <description>Used to specify integration directions output specifiers.</description>
   <encodeFunction>yes</encodeFunction>
   <visibility>public</visibility>
   <validator>yes</validator>
   <entry label="radial"                    />
   <entry label="lineOfSight"               />
   <entry label="lineOfSightInteriorAverage"/>
   <entry label="lambdaR"                   />
  </enumeration>
  !!]

contains

  subroutine Galactic_Structure_Radii_Definition_Decode(descriptors,specifiers,diskRequired,spheroidRequired,radiusVirialRequired,radiusScaleRequired)
    !!{
    Decode a set of radii descriptors and return the corresponding specifiers.
    !!}
    use :: Galactic_Structure_Options    , only : enumerationComponentTypeEncode   , enumerationMassTypeEncode, weightByLuminosity      , weightByMass , &
          &                                       weightIndexNull
    use :: Galacticus_Error              , only : Galacticus_Component_List        , Galacticus_Error_Report
    use :: Galacticus_Nodes              , only : defaultDarkMatterProfileComponent, defaultDiskComponent     , defaultSpheroidComponent, treeNode
    use :: ISO_Varying_String            , only : char                             , extract                  , operator(==)            , assignment(=), &
         &                                        operator(//)
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    use :: String_Handling               , only : String_Count_Words               , String_Split_Words       , char
    implicit none
    type     (varying_string ), intent(in   ), dimension(:)              :: descriptors
    type     (radiusSpecifier), intent(inout), dimension(:), allocatable :: specifiers
    logical                   , intent(  out)                            :: diskRequired        , spheroidRequired   , &
         &                                                                  radiusVirialRequired, radiusScaleRequired
    type     (varying_string  )              , dimension(5)              :: radiusDefinition
    type     (varying_string  )              , dimension(3)              :: fractionDefinition
    type     (varying_string  )              , dimension(2)              :: weightingDefinition
    type     (varying_string  )                                          :: valueDefinition     , message
    character(len=20          )                                          :: fractionLabel       , radiusLabel
    integer                                                              :: i                   , radiiCount         , &
         &                                                                  countComponents

    diskRequired        =.false.
    spheroidRequired    =.false.
    radiusVirialRequired=.false.
    radiusScaleRequired =.false.
    radiiCount          =size(descriptors)
    allocate(specifiers(radiiCount))
    do i=1,radiiCount
       specifiers(i)%name=descriptors(i)
       countComponents=String_Count_Words(char(descriptors(i)),':',bracketing="{}")
       call String_Split_Words(radiusDefinition,char(descriptors(i)),':',bracketing="{}")
       ! Detect cases which specify radius via a mass or light fraction. In either case, extract the fraction.
       valueDefinition=radiusDefinition(1)
       if (extract(valueDefinition,1,21) == 'galacticLightFraction') then
          call String_Split_Words(fractionDefinition,char(valueDefinition),'{}')
          radiusDefinition(1)='galacticLightFraction'
       end if
       if (extract(valueDefinition,1,20) == 'galacticMassFraction' ) then
          call String_Split_Words(fractionDefinition,char(valueDefinition),'{}')
          radiusDefinition(1)='galacticMassFraction'
       end if
       ! Parse the radius definition.
       select case (char(radiusDefinition(1)))
       case ('radius'                )
          specifiers(i)%type=radiusTypeRadius
       case ('virialRadius'          )
          specifiers(i)%type=radiusTypeVirialRadius
          radiusVirialRequired         =.true.
       case ('darkMatterScaleRadius' )
          specifiers(i)%type=radiusTypeDarkMatterScaleRadius
          radiusScaleRequired=.true.
          if (.not.defaultDarkMatterProfileComponent%scaleIsGettable         ())                                                   &
               & call Galacticus_Error_Report                                                                                      &
               &      (                                                                                                            &
               &       'dark matter profile scale radius is not gettable.'//                                                       &
               &        Galacticus_Component_List(                                                                                 &
               &                                  'darkMatterProfile'                                                           ,  &
               &                                   defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)   &
               &                                  )                                                                             // &
               &       {introspection:location}                                                                                    &
               &      )
       case ('diskRadius'            )
          specifiers(i)%type=radiusTypeDiskRadius
          diskRequired                 =.true.
          if (.not.defaultDiskComponent             %radiusIsGettable        ())                                                   &
               & call Galacticus_Error_Report                                                                                      &
               &(                                                                                                                  &
               &                              'disk radius is not gettable.'//                                                     &
               &        Galacticus_Component_List(                                                                                 &
               &                                  'disk'                                                                        ,  &
               &                                   defaultDiskComponent    %        radiusAttributeMatch(requireGettable=.true.)   &
               &                                  )                                                                             // &
               &       {introspection:location}                                                                                    &
               &                             )
       case ('spheroidRadius'        )
          specifiers(i)%type=radiusTypeSpheroidRadius
          spheroidRequired             =.true.
          if (.not.defaultSpheroidComponent         %radiusIsGettable        ())                                                   &
               & call Galacticus_Error_Report                                                                                      &
               &(                                                                                                                  &
               &                              'spheroid radius is not gettable.'//                                                 &
               &        Galacticus_Component_List(                                                                                 &
               &                                  'spheroid'                                                                    ,  &
               &                                   defaultSpheroidComponent%        radiusAttributeMatch(requireGettable=.true.)   &
               &                                  )                                                                             // &
               &       {introspection:location}                                                                                    &
               &                             )
       case ('diskHalfMassRadius'    )
          specifiers(i)%type=radiusTypeDiskHalfMassRadius
          diskRequired                 =.true.
          if (.not.defaultDiskComponent             %halfMassRadiusIsGettable())                                                   &
               & call Galacticus_Error_Report                                                                                      &
               &(                                                                                                                  &
               &                              'disk half-mass radius is not gettable.'//                                           &
               &        Galacticus_Component_List(                                                                                 &
               &                                  'disk'                                                                        ,  &
               &                                   defaultDiskComponent    %halfMassRadiusAttributeMatch(requireGettable=.true.)   &
               &                                  )                                                                             // &
               &       {introspection:location}                                                                                    &
               &                             )
       case ('spheroidHalfMassRadius')
          specifiers(i)%type=radiusTypeSpheroidHalfMassRadius
          spheroidRequired   =.true.
          if (.not.defaultSpheroidComponent         %halfMassRadiusIsGettable())                                                   &
               & call Galacticus_Error_Report                                                                                      &
               &(                                                                                                                  &
               &                              'spheroid half-mass radius is not gettable.'//                                       &
               &        Galacticus_Component_List(                                                                                 &
               &                                  'spheroid'                                                                    ,  &
               &                                   defaultSpheroidComponent%halfMassRadiusAttributeMatch(requireGettable=.true.)   &
               &                                  )                                                                             // &
               &       {introspection:location}                                                                                    &
               &                             )
       case ('galacticMassFraction'  )
          specifiers(i)%type=radiusTypeGalacticMassFraction
          fractionLabel=fractionDefinition(2)
          read (fractionLabel,*) specifiers(i)%fraction
          specifiers(i)%weightBy     =weightByMass
          specifiers(i)%weightByIndex=weightIndexNull
       case ('galacticLightFraction' )
          specifiers(i)%type=radiusTypeGalacticLightFraction
          fractionLabel=fractionDefinition(2)
          read (fractionLabel,*) specifiers(i)%fraction
          specifiers(i)%weightBy      =weightByLuminosity
          specifiers(i)%weightByIndex=unitStellarLuminosities%index(fractionDefinition(3))
       case default
          call Galacticus_Error_Report('unrecognized radius specifier "'//char(radiusDefinition(1))//'"'//{introspection:location})
       end select
       specifiers(i)%component=enumerationComponentTypeEncode(char(radiusDefinition(2)),includesPrefix=.false.)
       specifiers(i)%mass     =enumerationMassTypeEncode     (char(radiusDefinition(3)),includesPrefix=.false.)
       ! Detect cases which specify the weighting for integrals over the velocity dispersion.
       if (countComponents == 5) then
          valueDefinition=radiusDefinition(4)
          if     (                                                &
               &   extract(valueDefinition,1,11) == 'lineOfSight' &
               &  .or.                                            &
               &   extract(valueDefinition,1, 7) == 'lambdaR'     &
               & ) then
             call String_Split_Words(weightingDefinition,char(valueDefinition),'{}')
             radiusDefinition(4)=weightingDefinition(1)
             if (weightingDefinition(2) == "mass" .or. weightingDefinition(2) == "") then
                specifiers(i)%integralWeightBy     =weightByMass
                specifiers(i)%integralWeightByIndex=weightIndexNull
             else
                specifiers(i)%integralWeightBy     =weightByLuminosity
                specifiers(i)%integralWeightByIndex=unitStellarLuminosities%index(weightingDefinition(2))
             end if
          end if
          ! Parse the direction definition.
          select case (char(radiusDefinition(4)))
          case ('radial'                    )
             specifiers(i)%direction=directionRadial
          case ('lineOfSight'               )
             specifiers(i)%direction=directionLineOfSight
          case ('lineOfSightInteriorAverage')
             specifiers(i)%direction=directionLineOfSightInteriorAverage
          case ('lambdaR'                   )
             specifiers(i)%direction=directionLambdaR
          case default
             message='unrecognized direction specifier: "'//radiusDefinition(4)//'"'
             message=message//char(10)//'available specifiers are:'
             message=message//char(10)//' --> radial'
             message=message//char(10)//' --> lineOfSight'
             message=message//char(10)//' --> lineOfSightInteriorAverage'
             message=message//char(10)//' --> lambdaR'
             call Galacticus_Error_Report(message//{introspection:location})
          end select
       end if
       ! Get the numerical radius.
       radiusLabel=radiusDefinition(countComponents)
       read (radiusLabel,*) specifiers(i)%value
    end do
    return
  end subroutine Galactic_Structure_Radii_Definition_Decode

end module Galactic_Structure_Radii_Definitions
