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

!!{RST
Contains a module which provides parsing and evaluation of radii definitions used in output specifiers.
!!}

module Galactic_Structure_Radii_Definitions
  !!{RST
  Provides parsing and evaluation of radii definitions used in output specifiers.
  !!}
  use :: ISO_Varying_String        , only : varying_string
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType  , enumerationMassTypeType, enumerationWeightByType
  use :: Galacticus_Nodes          , only : nodeComponentDarkMatterProfile, nodeComponentDisk      , nodeComponentHotHalo   , &
       &                                    nodeComponentNSC              , nodeComponentSatellite , nodeComponentSpheroid  , &
       &                                    treeNode
  private
  public :: radiusSpecifier, radiusDefinitions, radiusResolver

  !!{RST
  Sentinel value returned in place of a radius which is undefined---for example, a radius specified as a fraction of a quantity
  which does not exist in the node in question. A large negative value is used since no physical radius can ever be negative, and
  since the code is compiled with ``-ffinite-math-only`` (so that ``NaN`` and ``∞`` are not usable as sentinels) and with
  ``-ffpe-trap=overflow`` (so that consumers **must** test for this sentinel *before* performing any arithmetic on the radius).
  !!}
  double precision, parameter, public :: radiusUndefined=-huge(0.0d0)

  !![
  <enumeration docformat="rst">
   <name>radiusType</name>
   <description>
   Used to specify radii types used in output specifiers.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="radius"                           description="Radii are specified absolutely, in units of Mpc"                                            />
   <entry label="virialRadius"                     description="Radii are specified in units of the virial radius"                                          />
   <entry label="hotHaloOuterRadius"               description="Radii are specified in units of the outer radius of the hot halo"                           />
   <entry label="darkMatterScaleRadius"            description="Radii are specified in units of the dark matter profile scale radius"                       />
   <entry label="diskRadius"                       description="Radii are specified in units of the disk scale radius"                                      />
   <entry label="nuclearStarClusterRadius"         description="Radii are specified in units of the \gls{nsc} scale radius"                                 />
   <entry label="spheroidRadius"                   description="Radii are specified in units of the spheroid scale radius"                                  />
   <entry label="diskHalfMassRadius"               description="Radii are specified in units of the disk half-mass radius"                                  />
   <entry label="nuclearStarClusterHalfMassRadius" description="Radii are specified in units of the \gls{nsc} half-mass radius"                             />
   <entry label="spheroidHalfMassRadius"           description="Radii are specified in units of the spheroid half-mass radius"                              />
   <entry label="satelliteBoundMassFraction"       description="Radii are specified in units of the radius enclosing a fraction of the satellite bound mass"/>
   <entry label="galacticMassFraction"             description="Radii are specified in units of the radius enclosing a fraction of the galactic mass"       />
   <entry label="galacticLightFraction"            description="Radii are specified in units of the radius enclosing a fraction of the galactic light"      />
   <entry label="stellarMassFraction"              description="Radii are specified in units of the radius enclosing a fraction of the stellar mass"        />
  </enumeration>
  !!]

  !![
  <enumeration docformat="rst">
   <name>direction</name>
   <description>
   Used to specify the type of velocity dispersion in output specifiers.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="radial"                     description="The radial velocity dispersion is computed"                                                                  />
   <entry label="lineOfSight"                description="The line-of-sight velocity dispersion (mass- or light-weighted) is computed"                                 />
   <entry label="lineOfSightInteriorAverage" description="The line-of-sight velocity dispersion, averaged over all interior radii (mass- or light-weighted)"           />
   <entry label="lambdaR"                    description="The λᵣ parameter of Ensellem et al. (2007; https://ui.adsabs.harvard.edu/abs/2007MNRAS.379..401E) is computed"/>
  </enumeration>
  !!]

  type radiusSpecifier
     !!{RST
     Type used for specifying radii definitions used in output specifiers.
     !!}
     type            (varying_string              ) :: name
     type            (enumerationComponentTypeType) :: component
     type            (enumerationMassTypeType     ) :: mass
     type            (enumerationWeightByType     ) :: weightBy             , integralWeightBy
     type            (enumerationDirectionType    ) :: direction
     type            (enumerationRadiusTypeType   ) :: type
     integer                                        :: integralWeightByIndex, weightByIndex
     double precision                               :: fraction             , value
  end type radiusSpecifier

  type radiusDefinitions
     !!{RST
     Type used to hold a decoded set of radius definitions, together with a record of which node components are required in order
     to evaluate them. Radii are evaluated in a given node by first constructing a :galacticus-type:`radiusResolver` from this
     object.
     !!}
     type   (radiusSpecifier), allocatable, dimension(:) :: specifiers
     integer                                             :: count                     =0
     logical                                             :: hotHaloRequired           =.false., diskRequired        =.false., &
          &                                                 spheroidRequired          =.false., satelliteRequired   =.false., &
          &                                                 nuclearStarClusterRequired=.false., radiusVirialRequired=.false., &
          &                                                 radiusScaleRequired       =.false.
   contains
     procedure :: decode => radiusDefinitionsDecode
  end type radiusDefinitions

  type radiusResolver
     !!{RST
     Type used to evaluate a set of :galacticus-type:`radiusDefinitions` in a single :galacticus-type:`treeNode`. Constructing a
     resolver gathers, once, all of the node components and halo scales needed by the definitions; each radius is then evaluated by
     a call to the ``evaluate`` method. A resolver holds pointers into the node it was built for, so it must be constructed as a
     local variable of the routine which uses it (and never cached in a shared object), and must not outlive that node.
     !!}
     type            (radiusDefinitions             ), pointer :: definitions_       => null()
     type            (treeNode                      ), pointer :: node_              => null()
     class           (nodeComponentHotHalo          ), pointer :: hotHalo            => null()
     class           (nodeComponentDisk             ), pointer :: disk               => null()
     class           (nodeComponentSpheroid         ), pointer :: spheroid           => null()
     class           (nodeComponentNSC              ), pointer :: nuclearStarCluster => null()
     class           (nodeComponentSatellite        ), pointer :: satellite          => null()
     class           (nodeComponentDarkMatterProfile), pointer :: darkMatterProfile  => null()
     double precision                                          :: radiusVirial       =  0.0d0, fractionDarkMatter=1.0d0
   contains
     procedure :: evaluate => radiusResolverEvaluate
  end type radiusResolver

  interface radiusResolver
     !!{RST
     Constructors for the :galacticus-type:`radiusResolver` type.
     !!}
     module procedure radiusResolverConstructor
  end interface radiusResolver

contains

  subroutine radiusDefinitionsDecode(self,descriptors)
    !!{RST
    Decode a set of radii descriptors and store the corresponding specifiers.
    !!}
    use :: Galactic_Structure_Options    , only : enumerationComponentTypeEncode   , enumerationMassTypeEncode  , weightByLuminosity      , weightByMass       , &
          &                                       enumerationComponentTypeDescribe , enumerationMassTypeDescribe, weightIndexNull
    use :: Error                         , only : Component_List                   , Error_Report               , errorStatusSuccess
    use :: Galacticus_Nodes              , only : defaultDarkMatterProfileComponent, defaultDiskComponent       , defaultSpheroidComponent, defaultNSCComponent, &
          &                                       defaultHotHaloComponent
    use :: ISO_Varying_String            , only : char                             , extract                    , operator(==)            , assignment(=)      , &
          &                                       operator(//)
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    use :: String_Handling               , only : String_Count_Words               , String_Split_Words         , char
    implicit none
    class    (radiusDefinitions), intent(inout), target       :: self
    type     (varying_string   ), intent(in   ), dimension(:) :: descriptors
    type     (radiusSpecifier  ), pointer      , dimension(:) :: specifiers
    type     (varying_string   )               , dimension(5) :: radiusDefinition
    type     (varying_string   )               , dimension(3) :: fractionDefinition
    type     (varying_string   )               , dimension(2) :: weightingDefinition
    type     (varying_string   )                              :: valueDefinition     , message
    character(len=20           )                              :: fractionLabel       , radiusLabel
    character(len=27           )                              :: extract_
    integer                                                   :: i                   , radiiCount         , &
         &                                                       countComponents     , status

    radiiCount=size(descriptors)
    self%count=radiiCount
    allocate(self%specifiers(radiiCount))
    ! Alias the specifiers array so that the parsing logic below need not repeat the `self%` prefix throughout.
    specifiers => self%specifiers
    do i=1,radiiCount
       specifiers(i)%name=descriptors(i)
       countComponents=String_Count_Words(char(descriptors(i)),':',bracketing="{}")
       if     (                     &
            &   countComponents < 4 &
            &  .or.                 &
            &   countComponents > 5 &
            & )  then
          message='radius specifier must have 4 (or, optionally, 5 if a direction is permitted) elements separated by `:`'
          call reportSpecifierError(specifiers(i)%name,message)
       end if
       call String_Split_Words(radiusDefinition,char(descriptors(i)),':',bracketing="{}")
       ! Default weights.
       specifiers(i)%weightBy             =weightByMass
       specifiers(i)%weightByIndex        =weightIndexNull
       specifiers(i)%integralWeightBy     =weightByMass
       specifiers(i)%integralWeightByIndex=weightIndexNull
       ! Detect cases which specify radius via a mass or light fraction. In either case, extract the fraction.
       extract_=extract(radiusDefinition(1),1,27)
       if (extract_(1:22) == 'galacticLightFraction{'     ) call extractFraction(specifiers(i)%name,radiusDefinition(1),22,fractionDefinition)
       if (extract_(1:21) == 'galacticMassFraction{'      ) call extractFraction(specifiers(i)%name,radiusDefinition(1),21,fractionDefinition)
       if (extract_(1:27) == 'satelliteBoundMassFraction{') call extractFraction(specifiers(i)%name,radiusDefinition(1),27,fractionDefinition)
       if (extract_(1:20) == 'stellarMassFraction{'       ) call extractFraction(specifiers(i)%name,radiusDefinition(1),20,fractionDefinition)
       ! Parse the radius definition.
       select case (char(radiusDefinition(1)))
       case ('radius'                          )
          specifiers(i)%type=radiusTypeRadius
       case ('virialRadius'                    )
          specifiers(i)%type=radiusTypeVirialRadius
          self%radiusVirialRequired      =.true.
       case ('hotHaloOuterRadius'              )
          specifiers(i)%type=radiusTypeHotHaloOuterRadius
          self%hotHaloRequired           =.true.
          if (.not.defaultHotHaloComponent          %outerRadiusIsGettable   ())                                        &
               & call Error_Report                                                                                      &
               &(                                                                                                       &
               &                              'hot halo outer radius is not gettable.'//                                &
               &        Component_List(                                                                                 &
               &                       'hotHalo'                                                                     ,  &
               &                        defaultHotHaloComponent %   outerRadiusAttributeMatch(requireGettable=.true.)   &
               &                       )                                                                             // &
               &       {introspection:location}                                                                         &
               &                             )
       case ('darkMatterScaleRadius'           )
          specifiers(i)%type=radiusTypeDarkMatterScaleRadius
          self%radiusScaleRequired       =.true.
          if (.not.defaultDarkMatterProfileComponent%scaleIsGettable         ())                                        &
               & call Error_Report                                                                                      &
               &      (                                                                                                 &
               &       'dark matter profile scale radius is not gettable.'//                                            &
               &        Component_List(                                                                                 &
               &                       'darkMatterProfile'                                                           ,  &
               &                        defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)   &
               &                       )                                                                             // &
               &       {introspection:location}                                                                         &
               &      )
       case ('diskRadius'                      )
          specifiers(i)%type=radiusTypeDiskRadius
          self%diskRequired              =.true.
          if (.not.defaultDiskComponent             %radiusIsGettable        ())                                        &
               & call Error_Report                                                                                      &
               &(                                                                                                       &
               &                              'disk radius is not gettable.'//                                          &
               &        Component_List(                                                                                 &
               &                       'disk'                                                                        ,  &
               &                        defaultDiskComponent    %        radiusAttributeMatch(requireGettable=.true.)   &
               &                       )                                                                             // &
               &       {introspection:location}                                                                         &
               &                             )
       case ('spheroidRadius'                  )
          specifiers(i)%type=radiusTypeSpheroidRadius
          self%spheroidRequired          =.true.
          if (.not.defaultSpheroidComponent         %radiusIsGettable        ())                                        &
               & call Error_Report                                                                                      &
               &(                                                                                                       &
               &                              'spheroid radius is not gettable.'//                                      &
               &        Component_List(                                                                                 &
               &                       'spheroid'                                                                    ,  &
               &                        defaultSpheroidComponent%        radiusAttributeMatch(requireGettable=.true.)   &
               &                       )                                                                             // &
               &       {introspection:location}                                                                         &
               &                             )
       case ('nuclearStarClusterRadius'        )
          specifiers(i)%type=radiusTypeNuclearStarClusterRadius
          self%nuclearStarClusterRequired=.true.
          if (.not.defaultNSCComponent             %radiusIsGettable        ())                                         &
               & call Error_Report                                                                                      &
               &(                                                                                                       &
               &                              'nuclear star cluster radius is not gettable.'//                          &
               &        Component_List(                                                                                 &
               &                       'NSC'                                                                        ,   &
               &                        defaultNSCComponent    %        radiusAttributeMatch(requireGettable=.true.)    &
               &                       )                                                                             // &
               &       {introspection:location}                                                                         &
               &                             )

       case ('diskHalfMassRadius'              )
          specifiers(i)%type=radiusTypeDiskHalfMassRadius
          self%diskRequired              =.true.
          if (.not.defaultDiskComponent             %halfMassRadiusIsGettable())                                        &
               & call Error_Report                                                                                      &
               &(                                                                                                       &
               &                              'disk half-mass radius is not gettable.'//                                &
               &        Component_List(                                                                                 &
               &                       'disk'                                                                        ,  &
               &                        defaultDiskComponent    %halfMassRadiusAttributeMatch(requireGettable=.true.)   &
               &                       )                                                                             // &
               &       {introspection:location}                                                                         &
               &                             )
       case ('spheroidHalfMassRadius'          )
          specifiers(i)%type=radiusTypeSpheroidHalfMassRadius
          self%spheroidRequired          =.true.
          if (.not.defaultSpheroidComponent         %halfMassRadiusIsGettable())                                        &
               & call Error_Report                                                                                      &
               &(                                                                                                       &
               &                              'spheroid half-mass radius is not gettable.'//                            &
               &        Component_List(                                                                                 &
               &                       'spheroid'                                                                    ,  &
               &                        defaultSpheroidComponent%halfMassRadiusAttributeMatch(requireGettable=.true.)   &
               &                       )                                                                             // &
               &       {introspection:location}                                                                         &
               &                             )
       case ('nuclearStarClusterHalfMassRadius')
          specifiers(i)%type=radiusTypeNuclearStarClusterHalfMassRadius
          self%nuclearStarClusterRequired=.true.
          if (.not.defaultNSCComponent         %halfMassRadiusIsGettable())                                             &
               & call Error_Report                                                                                      &
               &(                                                                                                       &
               &                              'nuclear star cluster half-mass radius is not gettable.'//                &
               &        Component_List(                                                                                 &
               &                       'NSC'                                                                         ,  &
               &                        defaultNSCComponent     %halfMassRadiusAttributeMatch(requireGettable=.true.)   &
               &                       )                                                                             // &
               &       {introspection:location}                                                                         &
               &                             )

       case ('satelliteBoundMassFraction'      )
          specifiers(i)%type=radiusTypeSatelliteBoundMassFraction
          self%satelliteRequired         =.true.
          fractionLabel=fractionDefinition(2)
          read (fractionLabel,*,iostat=status) specifiers(i)%fraction
          if (status /= 0) then
             message='unable to parse numerical fraction'
             call reportSpecifierError(specifiers(i)%name,message,highlight=1,bracketed=.true.)
          end if
          specifiers(i)%weightBy     =weightByMass
          specifiers(i)%weightByIndex=weightIndexNull
       case ('galacticMassFraction'            )
          specifiers(i)%type=radiusTypeGalacticMassFraction
          fractionLabel=fractionDefinition(2)
          read (fractionLabel,*,iostat=status) specifiers(i)%fraction
          if (status /= 0) then
             message='unable to parse numerical fraction'
             call reportSpecifierError(specifiers(i)%name,message,highlight=1,bracketed=.true.)
          end if
          specifiers(i)%weightBy     =weightByMass
          specifiers(i)%weightByIndex=weightIndexNull
       case ('galacticLightFraction'           )
          specifiers(i)%type=radiusTypeGalacticLightFraction
          fractionLabel=fractionDefinition(2)
          read (fractionLabel,*,iostat=status) specifiers(i)%fraction
          if (status /= 0) then
             message='unable to parse numerical fraction'
             call reportSpecifierError(specifiers(i)%name,message,highlight=1,bracketed=.true.)
          end if
          specifiers(i)%weightBy      =weightByLuminosity
          specifiers(i)%weightByIndex=unitStellarLuminosities%index(fractionDefinition(3))
       case ('stellarMassFraction'             )
          specifiers(i)%type=radiusTypeStellarMassFraction
          fractionLabel=fractionDefinition(2)
          read (fractionLabel,*,iostat=status) specifiers(i)%fraction
          if (status /= 0) then
             message='unable to parse numerical fraction'
             call reportSpecifierError(specifiers(i)%name,message,highlight=1,bracketed=.true.)
          end if
          specifiers(i)%weightBy     =weightByMass
          specifiers(i)%weightByIndex=weightIndexNull
       case default
          message="unrecognized radius type"//char(10)//enumerationRadiusTypeDescribe()
          call reportSpecifierError(specifiers(i)%name,message,highlight=1,bracketed=.false.)
       end select
       specifiers(i)%component=enumerationComponentTypeEncode(char(radiusDefinition(2)),includesPrefix=.false.,status=status)
       if (status /= errorStatusSuccess) then
          message="unrecognized component type"//char(10)//enumerationComponentTypeDescribe()
          call reportSpecifierError(specifiers(i)%name,message,highlight=2,bracketed=.false.)
       end if
       specifiers(i)%mass     =enumerationMassTypeEncode     (char(radiusDefinition(3)),includesPrefix=.false.,status=status)
       if (status /= errorStatusSuccess) then
          message="unrecognized mass type"    //char(10)//enumerationMassTypeDescribe      ()
          call reportSpecifierError(specifiers(i)%name,message,highlight=3,bracketed=.false.)
       end if
       ! Detect cases which specify the weighting for integrals over the velocity dispersion.
       if (countComponents == 5) then
          valueDefinition=radiusDefinition(4)
          extract_       =extract(valueDefinition,1,11)
          if     (                                 &
               &   extract_(1:11) == 'lineOfSight' &
               &  .or.                             &
               &   extract_(1: 7) == 'lambdaR'     &
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
             message='unrecognized direction specifier'//char(10)//enumerationDirectionDescribe()
             call reportSpecifierError(specifiers(i)%name,message,highlight=4,bracketed=.false.)
          end select
       end if
       ! Get the numerical radius.
       radiusLabel=radiusDefinition(countComponents)
       read (radiusLabel,*,iostat=status) specifiers(i)%value
       if (status /= 0) then
          message='unable to parse numerical radius'
          call reportSpecifierError(specifiers(i)%name,message,highlight=countComponents,bracketed=.true.)
       end if
    end do
    return
  end subroutine radiusDefinitionsDecode

  function radiusResolverConstructor(definitions_,node,darkMatterHaloScale_,fractionDarkMatter,radiusVirialRequired) result(self)
    !!{RST
    Constructor for the :galacticus-type:`radiusResolver` type. Gathers, once, all of the node components and halo scales which
    are required in order to evaluate the given set of radius definitions in the given node.

    The optional {\normalfont \ttfamily fractionDarkMatter} argument gives the fraction of the bound mass of a satellite which is
    dark matter---it is required only if radii are specified in units of the satellite bound mass fraction. The optional
    {\normalfont \ttfamily radiusVirialRequired} argument forces the virial radius to be computed even if no radius definition
    requires it---it is used by consumers which need the virial radius for their own purposes.
    !!}
    implicit none
    type            (radiusResolver          )                          :: self
    type            (radiusDefinitions       ), intent(in   ), target    :: definitions_
    type            (treeNode                ), intent(inout), target    :: node
    class           (darkMatterHaloScaleClass), intent(inout)            :: darkMatterHaloScale_
    double precision                          , intent(in   ), optional  :: fractionDarkMatter
    logical                                   , intent(in   ), optional  :: radiusVirialRequired
    !![
    <optionalArgument name="fractionDarkMatter"   defaultsTo="1.0d0"          />
    <optionalArgument name="radiusVirialRequired" defaultsTo=".false."        />
    !!]

    self%definitions_      => definitions_
    self%node_             => node
    self%fractionDarkMatter=  fractionDarkMatter_
    if (definitions_%radiusVirialRequired.or.radiusVirialRequired_)                                     &
         & self%radiusVirial       =  darkMatterHaloScale_%radiusVirial(node                    )
    if (definitions_%hotHaloRequired              ) self%hotHalo            => node%hotHalo          ()
    if (definitions_%diskRequired                 ) self%disk               => node%disk             ()
    if (definitions_%spheroidRequired             ) self%spheroid           => node%spheroid         ()
    if (definitions_%nuclearStarClusterRequired   ) self%nuclearStarCluster => node%NSC              ()
    if (definitions_%satelliteRequired            ) self%satellite          => node%satellite        ()
    if (definitions_%radiusScaleRequired          ) self%darkMatterProfile  => node%darkMatterProfile()
    return
  end function radiusResolverConstructor

  subroutine radiusResolverEvaluate(self,indexRadius,radius,radiusScale)
    !!{RST
    Evaluate the radius corresponding to the {\normalfont \ttfamily indexRadius}\ :sup:`th` radius definition in the node to which
    this resolver is attached.

    The optional {\normalfont \ttfamily radiusScale} argument returns the scale radius on which the radius is based---that is, the
    quantity of which the specified radius is a multiple (for radii specified absolutely, in Mpc, the radius itself is
    returned). This is useful for consumers which must choose an outer radius for an integral.

    If the radius is undefined, {\normalfont \ttfamily radiusUndefined} is returned in both arguments. Consumers **must** test for
    this sentinel *before* performing any arithmetic on the returned radius, since the sentinel is a large negative number which
    will overflow (and so trap) under multiplication.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure_Options, only : componentTypeAll     , massTypeDark, massTypeGalactic, &
         &                                     massTypeStellar
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (radiusResolver       ), intent(inout)           :: self
    integer                                , intent(in   )           :: indexRadius
    double precision                       , intent(  out)           :: radius
    double precision                       , intent(  out), optional :: radiusScale
    class           (massDistributionClass), pointer                 :: massDistribution_
    type            (radiusSpecifier      ), pointer                 :: specifier_
    double precision                                                 :: radiusScale_     , mass

    specifier_ => self%definitions_%specifiers(indexRadius)
    radius     =  specifier_%value
    select case (specifier_%type%ID)
    case   (radiusTypeRadius                          %ID)
       radiusScale_=radius
    case   (radiusTypeVirialRadius                    %ID)
       radiusScale_=self%radiusVirial
    case   (radiusTypeDarkMatterScaleRadius           %ID)
       radiusScale_=self%darkMatterProfile %         scale()
    case   (radiusTypeHotHaloOuterRadius              %ID)
       radiusScale_=self%hotHalo           %   outerRadius()
    case   (radiusTypeDiskRadius                      %ID)
       radiusScale_=self%disk              %        radius()
    case   (radiusTypeSpheroidRadius                  %ID)
       radiusScale_=self%spheroid          %        radius()
    case   (radiusTypeNuclearStarClusterRadius        %ID)
       radiusScale_=self%nuclearStarCluster%        radius()
    case   (radiusTypeDiskHalfMassRadius              %ID)
       radiusScale_=self%disk              %halfMassRadius()
    case   (radiusTypeSpheroidHalfMassRadius          %ID)
       radiusScale_=self%spheroid          %halfMassRadius()
    case   (radiusTypeNuclearStarClusterHalfMassRadius%ID)
       radiusScale_=self%nuclearStarCluster%halfMassRadius()
    case   (radiusTypeSatelliteBoundMassFraction      %ID)
       mass              =  +self             %satellite%boundMass          (                                                &
            &                                                               )                                                &
            &               *self             %fractionDarkMatter
       massDistribution_ =>  self             %node_    %massDistribution   (                                                &
            &                                                                massType      =              massTypeDark    ,  &
            &                                                                componentType =              componentTypeAll,  &
            &                                                                weightBy      =specifier_%weightBy           ,  &
            &                                                                weightIndex   =specifier_%weightByIndex         &
            &                                                               )
       radiusScale_      =  +massDistribution_%radiusEnclosingMass(                                                          &
            &                                                      mass          =              mass                        &
            &                                                     )
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    case   (radiusTypeGalacticMassFraction            %ID)
       massDistribution_ =>  self             %node_    %massDistribution   (                                                &
            &                                                                massType      =              massTypeGalactic,  &
            &                                                                componentType =              componentTypeAll,  &
            &                                                                weightBy      =specifier_%weightBy           ,  &
            &                                                                weightIndex   =specifier_%weightByIndex         &
            &                                                               )
       radiusScale_      =  +massDistribution_%radiusEnclosingMass(                                                          &
            &                                                      massFractional=specifier_%fraction                       &
            &                                                     )
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    case   (radiusTypeGalacticLightFraction           %ID,  &
         &  radiusTypeStellarMassFraction             %ID)
       massDistribution_ =>  self             %node_    %massDistribution   (                                                &
            &                                                                massType      =              massTypeStellar ,  &
            &                                                                componentType =              componentTypeAll,  &
            &                                                                weightBy      =specifier_%weightBy           ,  &
            &                                                                weightIndex   =specifier_%weightByIndex         &
            &                                                               )
       radiusScale_      =  +massDistribution_%radiusEnclosingMass(                                                          &
            &                                                      massFractional=specifier_%fraction                       &
            &                                                     )
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    case default
       radiusScale_=radiusUndefined
       call Error_Report('unrecognized radius type'//{introspection:location})
    end select
    ! For radii specified absolutely there is no scale radius to multiply by - the radius is simply the value given.
    if (specifier_%type%ID /= radiusTypeRadius%ID) radius=radius*radiusScale_
    if (present(radiusScale)) radiusScale=radiusScale_
    return
  end subroutine radiusResolverEvaluate

  subroutine extractFraction(specifier,radiusDefinition,openAt,fractionDefinition)
    !!{RST
    Parse a fractional radius definition.
    !!}
    use :: ISO_Varying_String, only : extract           , char, len, index, assignment(=)
    use :: String_Handling   , only : String_Split_Words
    implicit none
    type   (varying_string), intent(in   )               :: specifier
    type   (varying_string), intent(inout)               :: radiusDefinition
    integer                , intent(in   )               :: openAt
    type   (varying_string), intent(  out), dimension(3) :: fractionDefinition
    type   (varying_string)                              :: message

    if (index(radiusDefinition,'}') /= len(radiusDefinition)) then
       message='missing `}`'
       call reportSpecifierError(specifier,message,highlight=1)
    end if
    call String_Split_Words(fractionDefinition,char(radiusDefinition),'{}')
    radiusDefinition=extract(radiusDefinition,1,openAt-1)
    return
  end subroutine extractFraction
  
  subroutine reportSpecifierError(specifier,message,highlight,bracketed)
    !!{RST
    Report an error in parsing a radius specifier.
    !!}
    use :: Display           , only : displayGreen      , displayRed        , displayReset
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : operator(//)      , var_str           , char        , extract, &
         &                            index
    use :: String_Handling   , only : String_Count_Words, String_Split_Words, String_Join
    implicit none
    type   (varying_string), intent(in   )               :: specifier       , message
    integer                , intent(in   ), optional     :: highlight
    logical                , intent(in   ), optional     ::  bracketed
    type   (varying_string)               , dimension(5) :: radiusDefinition
    type   (varying_string)                              :: specifier_
    integer                                              :: countComponents , indexBracket
    !![
    <optionalArgument name="bracketed" defaultsTo=".false."/>
    !!]
    
    if (present(highlight)) then
       countComponents=String_Count_Words(char(specifier),':',bracketing="{}")
       call String_Split_Words(radiusDefinition,char(specifier),':',bracketing="{}")
       indexBracket=index(radiusDefinition(highlight),'{')
       if (indexBracket /= 0) then
          if (bracketed_) then
             radiusDefinition(highlight)=              extract(radiusDefinition(highlight),1,indexBracket-1)//displayRed  ()//extract(radiusDefinition(highlight),indexBracket)//displayReset()
          else
             radiusDefinition(highlight)=displayRed()//extract(radiusDefinition(highlight),1,indexBracket-1)//displayReset()//extract(radiusDefinition(highlight),indexBracket)
          end if
       else
          radiusDefinition   (highlight)=displayRed()//        radiusDefinition(highlight)                  //displayReset()
       end if
       specifier_=String_Join(radiusDefinition(1:countComponents),":")
    else
       specifier_=specifier
    end if
    call Error_Report(var_str('Failed to parse radius specifier:')//char(10)//char(10)//'   '//specifier_//char(10)//char(10)//message//char(10)//char(10)//displayGreen()//'HELP:'//displayReset()//' See https://galacticus.readthedocs.io/en/latest/manuals/physics/definitions.html#manual-sec-radiusspecifiers for an explanation of radius specifier syntax'//{introspection:location})    
    return
  end subroutine reportSpecifierError
  
end module Galactic_Structure_Radii_Definitions
