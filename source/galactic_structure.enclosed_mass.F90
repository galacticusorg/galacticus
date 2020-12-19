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

!% Contains a module which implements calculations of the mass enclosed within a specified radius.

module Galactic_Structure_Enclosed_Masses
  !% Implements calculations of the mass enclosed within a specified radius.
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private
  public :: Galactic_Structure_Enclosed_Mass, Galactic_Structure_Radius_Enclosing_Mass, Galactic_Structure_Radius_Enclosing_Density

  ! Variables used in root finding.
  integer                             :: componentTypeShared, massTypeShared, weightByShared, &
       &                                 weightIndexShared
  double precision                    :: massRoot           , radiusShared  , densityRoot
  type            (treeNode), pointer :: activeNode
  !$omp threadprivate(massRoot,densityRoot,radiusShared,massTypeShared,componentTypeShared,weightByShared,weightIndexShared,activeNode)

contains

  double precision function Galactic_Structure_Enclosed_Mass(node,radius,componentType,massType,weightBy,weightIndex)
    !% Solve for the mass within a given radius, or the total mass if no radius is specified. Assumes that galactic structure has
    !% already been computed.
    use :: Galactic_Structure_Options, only : radiusLarge
    use :: Galacticus_Nodes          , only : optimizeForEnclosedMassSummation, reductionSummation, treeNode
    !# <include directive="enclosedMassTask" type="moduleUse">
    include 'galactic_structure.enclosed_mass.tasks.modules.inc'
    !# </include>
    implicit none
    type            (treeNode               ), intent(inout)           :: node
    integer                                  , intent(in   ), optional :: componentType        , massType, weightBy, &
         &                                                                weightIndex
    double precision                         , intent(in   ), optional :: radius
    procedure       (Component_Enclosed_Mass), pointer                 :: componentEnclosedMass
    double precision                                                   :: componentMass

    ! Set default options.
    call Galactic_Structure_Enclosed_Mass_Defaults(componentType,massType,weightBy,weightIndex)
    ! Determine which radius to use.
    if (present(radius)) then
       radiusShared=radius
    else
       radiusShared=radiusLarge
    end if
    ! Compute the contribution from components directly, by mapping a function over all components.
    componentEnclosedMass => Component_Enclosed_Mass
    Galactic_Structure_Enclosed_Mass=node%mapDouble0(componentEnclosedMass,reductionSummation,optimizeFor=optimizeForEnclosedMassSummation)
    ! Call routines to supply the masses for all components.
    !# <include directive="enclosedMassTask" type="functionCall" functionType="function" returnParameter="componentMass">
    !#  <functionArgs>node,radiusShared,componentTypeShared,massTypeShared,weightByShared,weightIndexShared</functionArgs>
    !#  <onReturn>Galactic_Structure_Enclosed_Mass=Galactic_Structure_Enclosed_Mass+componentMass</onReturn>
    include 'galactic_structure.enclosed_mass.tasks.inc'
    !# </include>
    return
  end function Galactic_Structure_Enclosed_Mass

  double precision function Galactic_Structure_Radius_Enclosing_Mass(node,mass,fractionalMass,componentType,massType,weightBy,weightIndex)
    !% Return the radius enclosing a given mass (or fractional mass) in {\normalfont \ttfamily node}.
    use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale       , darkMatterHaloScaleClass
    use :: Dark_Matter_Profiles      , only : darkMatterProfile         , darkMatterProfileClass
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo     , massTypeDark
    use :: Galacticus_Display        , only : Galacticus_Display_Message, verbosityWarn
    use :: Galacticus_Error          , only : Galacticus_Error_Report
    use :: ISO_Varying_String        , only : varying_string            , assignment(=)                , operator(//)
    use :: Kind_Numbers              , only : kind_int8
    use :: Root_Finder               , only : rangeExpandMultiplicative , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    use :: String_Handling           , only : operator(//)
    implicit none
    type            (treeNode                ), intent(inout), target   :: node
    integer                                   , intent(in   ), optional :: componentType                        , massType   , &
         &                                                                 weightBy                             , weightIndex
    double precision                          , intent(in   ), optional :: fractionalMass                       , mass
    class           (darkMatterHaloScaleClass), pointer                 :: darkMatterHaloScale_
    class           (darkMatterProfileClass  ), pointer                 :: darkMatterProfile_
    type            (rootFinder              ), save                    :: finder
    !$omp threadprivate(finder)
    double precision                          , save                    :: radiusPrevious          =-huge(0.0d0)
    integer         (kind_int8               ), save                    :: uniqueIDPrevious        =-1_kind_int8
    !$omp threadprivate(radiusPrevious,uniqueIDPrevious)
    double precision                                                    :: radiusGuess
    type            (varying_string          )                          :: message
    character       (len=11                  )                          :: massLabel

     ! Set default options.
    call Galactic_Structure_Enclosed_Mass_Defaults(componentType,massType,weightBy,weightIndex)
    ! Determine what mass to use.
    if (present(mass)) then
       if (present(fractionalMass)) call Galacticus_Error_Report('only one mass or fractionalMass can be specified'//{introspection:location})
       massRoot=mass
    else if (present(fractionalMass)) then
       massRoot=fractionalMass*Galactic_Structure_Enclosed_Mass(node,componentType=componentTypeShared,massType=massTypeShared,weightBy=weightByShared,weightIndex=weightIndexShared)
    else
       call Galacticus_Error_Report('either mass or fractionalMass must be specified'//{introspection:location})
    end if
    if (massRoot <= 0.0d0) then
       Galactic_Structure_Radius_Enclosing_Mass=0.0d0
       return
    end if
    activeNode => node
    ! If dark matter component is queried and its density profile is unaffected by baryons, compute the radius from dark
    ! matter profile. Otherwise, find the radius numerically.
    if     (                                              &
         &   componentTypeShared == componentTypeDarkHalo &
         &  .or.                                          &
         &   massTypeShared      == massTypeDark          &
         & ) then
       darkMatterProfile_ => darkMatterProfile()
       Galactic_Structure_Radius_Enclosing_Mass=darkMatterProfile_%radiusEnclosingMass(activeNode,massRoot)
    else
       ! Initialize our root finder.
       if (.not.finder%isInitialized()) then
          call finder%rangeExpand (                                                             &
               &                   rangeExpandDownward          =0.5d0                        , &
               &                   rangeExpandUpward            =2.0d0                        , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
          call finder%rootFunction(Enclosed_Mass_Root                              )
          call finder%tolerance   (toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
       end if
       ! Solve for the radius.
       if (Enclosed_Mass_Root(0.0d0) >= 0.0d0) then
          message='Enclosed mass in galaxy (ID='
          write (massLabel,'(e10.4)') Galactic_Structure_Enclosed_Mass(activeNode,0.0d0,componentTypeShared,massTypeShared,weightByShared,weightIndexShared)
          message=message//node%index()//') seems to be finite ('//trim(massLabel)
          write (massLabel,'(e10.4)') massRoot
          message=message//') at zero radius (was seeking '//trim(massLabel)
          message=message//') - returning zero radius.'
          call Galacticus_Display_Message(message,verbosityWarn)
          Galactic_Structure_Radius_Enclosing_Mass=0.0d0
          return
       end if
       if (node%uniqueID() == uniqueIDPrevious) then
          radiusGuess          =  radiusPrevious
       else
          darkMatterHaloScale_ => darkMatterHaloScale              (        )
          radiusGuess          =  darkMatterHaloScale_%virialRadius(node)
       end if
       Galactic_Structure_Radius_Enclosing_Mass=finder%find(rootGuess=radiusGuess)
       uniqueIDPrevious                        =node%uniqueID()
       radiusPrevious                          =Galactic_Structure_Radius_Enclosing_Mass
    end if
    return
  end function Galactic_Structure_Radius_Enclosing_Mass

  double precision function Galactic_Structure_Radius_Enclosing_Density(node,density,densityContrast,componentType,massType,weightBy,weightIndex)
    !% Return the radius enclosing a given density (or density contrast) in {\normalfont \ttfamily node}.
    use :: Cosmology_Functions    , only : cosmologyFunctions       , cosmologyFunctionsClass
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale      , darkMatterHaloScaleClass
    use :: Galacticus_Error       , only : Galacticus_Error_Report
    use :: Galacticus_Nodes       , only : nodeComponentBasic       , treeNode
    use :: Root_Finder            , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    type            (treeNode                ), intent(inout), target   :: node
    integer                                   , intent(in   ), optional :: componentType       , massType       , weightBy, weightIndex
    double precision                          , intent(in   ), optional :: density             , densityContrast
    class           (darkMatterHaloScaleClass)               , pointer  :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass )               , pointer  :: cosmologyFunctions_
    class           (nodeComponentBasic      )               , pointer  :: basic
    type            (rootFinder              ), save                    :: finder
    !$omp threadprivate(finder)

    ! Set default options.
    call Galactic_Structure_Enclosed_Mass_Defaults(componentType,massType,weightBy,weightIndex)
    ! Determine what mass to use.
    if (present(density)) then
       if (present(densityContrast)) call Galacticus_Error_Report('only one density or densityContrast can be specified'//{introspection:location})
       densityRoot=density
    else if (present(densityContrast)) then
       cosmologyFunctions_ => cosmologyFunctions()
       basic               => node%basic    ()
       densityRoot=densityContrast*cosmologyFunctions_%matterDensityEpochal(time=basic%time())
    else
       call Galacticus_Error_Report('either density or densityContrast must be specified'//{introspection:location})
    end if
    if (densityRoot <= 0.0d0) then
       Galactic_Structure_Radius_Enclosing_Density=0.0d0
       return
    end if
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rangeExpand (                                                             &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
       call finder%rootFunction(Enclosed_Density_Root                           )
       call finder%tolerance   (toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
    end if
    ! Solve for the radius.
    activeNode           => node
    darkMatterHaloScale_ => darkMatterHaloScale()
    massRoot             =  0.0d0
    Galactic_Structure_Radius_Enclosing_Density=finder%find(rootGuess=darkMatterHaloScale_%virialRadius(node))
    return
  end function Galactic_Structure_Radius_Enclosing_Density

  subroutine Galactic_Structure_Enclosed_Mass_Defaults(componentType,massType,weightBy,weightIndex)
    !% Set the default values for options in the enclosed mass functions.
    use :: Galactic_Structure_Options, only : componentTypeAll       , massTypeAll, weightByLuminosity, weightByMass
    use :: Galacticus_Error          , only : Galacticus_Error_Report
    implicit none
    integer, intent(in   ), optional :: componentType, massType, weightBy, weightIndex

    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeShared=massType
    else
       massTypeShared=massTypeAll
    end if
    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Determine which weighting to use.
    if (present(weightBy)) then
       weightByShared=weightBy
       select case (weightByShared)
       case (weightByLuminosity)
          if (.not.present(weightIndex)) call Galacticus_Error_Report('weightIndex should be specified for luminosity weighting'//{introspection:location})
          weightIndexShared=weightIndex
       end select
    else
       weightByShared=weightByMass
    end if
    return
  end subroutine Galactic_Structure_Enclosed_Mass_Defaults

  double precision function Enclosed_Mass_Root(radius)
    !% Root function used in solving for the radius that encloses a given mass.
    double precision, intent(in   ) :: radius

    ! Evaluate the root function.
    Enclosed_Mass_Root=Galactic_Structure_Enclosed_Mass(activeNode,radius,componentTypeShared,massTypeShared,weightByShared&
         &,weightIndexShared)-massRoot
    return
  end function Enclosed_Mass_Root

  double precision function Enclosed_Density_Root(radius)
    !% Root function used in solving for the radius that encloses a given density.
    use :: Numerical_Constants_Math, only : Pi
    double precision, intent(in   ) :: radius

    ! Evaluate the root function.
    Enclosed_Density_Root=3.0d0*Enclosed_Mass_Root(radius)/4.0d0/Pi/radius**3-densityRoot
  end function Enclosed_Density_Root

  double precision function Component_Enclosed_Mass(component)
    !% Unary function returning the enclosed mass in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    Component_Enclosed_Mass=component%enclosedMass(radiusShared,componentTypeShared,massTypeShared&
         &,weightByShared,weightIndexShared)
    return
  end function Component_Enclosed_Mass

end module Galactic_Structure_Enclosed_Masses
