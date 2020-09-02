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

!% Contains a module which implements a property extractor class for the mass and radii of spheres are specified density contrast.

  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Cosmology_Parameters   , only : cosmologyParameters, cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass
  use :: Galacticus_Nodes       , only : nodeComponentBasic , treeNode

  !# <nodePropertyExtractor name="nodePropertyExtractorDensityContrasts">
  !#  <description>A property extractor class for the mass and radii of spheres are specified density contrast.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorDensityContrasts
     !% A property extractor class for the mass and radii of spheres are specified density contrast..
     private
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_  => null()
     class           (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_ => null()
     integer                                                               :: elementCount_                 , countDensityContrasts, &
          &                                                                   massTypeSelected
     logical                                                               :: darkMatterOnly
     double precision                                                      :: densityReference
     double precision                          , allocatable, dimension(:) :: densityContrasts
   contains
     final     ::                 densityContrastsDestructor
     procedure :: elementCount => densityContrastsElementCount
     procedure :: extract      => densityContrastsExtract
     procedure :: names        => densityContrastsNames
     procedure :: descriptions => densityContrastsDescriptions
     procedure :: unitsInSI    => densityContrastsUnitsInSI
     procedure :: type         => densityContrastsType
  end type nodePropertyExtractorDensityContrasts

  interface nodePropertyExtractorDensityContrasts
     !% Constructors for the ``densityContrasts'' output analysis class.
     module procedure densityContrastsConstructorParameters
     module procedure densityContrastsConstructorInternal
  end interface nodePropertyExtractorDensityContrasts

  ! Module-scope variables used in root finding.
  type            (treeNode                             ), pointer :: densityContrastsNode
  class           (nodeComponentBasic                   ), pointer :: densityContrastsBasic
  class           (nodePropertyExtractorDensityContrasts), pointer :: densityContrastsSelf
  double precision                                                 :: densityContrastsMeanTarget
  !$omp threadprivate(densityContrastsNode,densityContrastsBasic,densityContrastsMeanTarget)

contains

  function densityContrastsConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily densityContrasts} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorDensityContrasts)                              :: self
    type            (inputParameters                      ), intent(inout)               :: parameters
    class           (cosmologyParametersClass             ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass              ), pointer                     :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass             ), pointer                     :: darkMatterHaloScale_
    double precision                                       , allocatable  , dimension(:) :: densityContrasts
    logical                                                                              :: darkMatterOnly

    allocate(densityContrasts(parameters%count('densityContrasts')))
    !# <inputParameter>
    !#   <name>densityContrasts</name>
    !#   <description>A list of density contrasts at which to output data.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>darkMatterOnly</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not density contrast data should be computed using the dark matter component alone.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=nodePropertyExtractorDensityContrasts(densityContrasts,darkMatterOnly,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function densityContrastsConstructorParameters

  function densityContrastsConstructorInternal(densityContrasts,darkMatterOnly,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily densityContrasts} property extractor class.
    use :: Galactic_Structure_Options, only : massTypeAll, massTypeDark
    implicit none
    type            (nodePropertyExtractorDensityContrasts)                              :: self
    class           (cosmologyFunctionsClass              ), intent(in   ), target       :: cosmologyFunctions_
    class           (cosmologyParametersClass             ), intent(in   ), target       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass             ), intent(in   ), target       :: darkMatterHaloScale_
    double precision                                       , intent(in   ), dimension(:) :: densityContrasts
    logical                                                , intent(in   )               :: darkMatterOnly
    !# <constructorAssign variables="densityContrasts, darkMatterOnly, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_"/>

    self%countDensityContrasts=  size(densityContrasts)
    self%elementCount_        =2*size(densityContrasts)
     select case (darkMatterOnly)
    case (.true.)
       self%massTypeSelected=massTypeDark
       self%densityReference=(self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon())*self%cosmologyParameters_%densityCritical()
    case (.false.)
       self%massTypeSelected=massTypeAll
       self%densityReference= self%cosmologyParameters_%OmegaMatter()                                         *self%cosmologyParameters_%densityCritical()
    end select
    return
  end function densityContrastsConstructorInternal

  subroutine densityContrastsDestructor(self)
    !% Destructor for the {\normalfont \ttfamily densityContrasts} property extractor class.
    implicit none
    type(nodePropertyExtractorDensityContrasts), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%cosmologyFunctions_" />
    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine densityContrastsDestructor

  integer function densityContrastsElementCount(self,time)
    !% Return the number of elements in the {\normalfont \ttfamily densityContrasts} property extractors.
    implicit none
    class           (nodePropertyExtractorDensityContrasts), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    densityContrastsElementCount=self%elementCount_
    return
  end function densityContrastsElementCount

  function densityContrastsExtract(self,node,time,instance)
    !% Implement a last isolated redshift output analysis.
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options        , only : componentTypeAll
    use :: Root_Finder                       , only : rangeExpandMultiplicative       , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    double precision                                        , dimension(:) , allocatable :: densityContrastsExtract
    class           (nodePropertyExtractorDensityContrasts ), intent(inout), target      :: self
    type            (treeNode                              ), intent(inout), target      :: node
    double precision                                        , intent(in   )              :: time
    type            (multiCounter                          ), intent(inout), optional    :: instance
    double precision                                        , parameter                  :: toleranceAbsolute      =0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder                            ), save                       :: finder
    !$omp threadprivate(finder)
    integer                                                                              :: i
    double precision                                                                     :: enclosedMass                 , radius
    !$GLC attributes unused :: time, instance

    allocate(densityContrastsExtract(self%elementCount_))
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rangeExpand (                                                             &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
       call finder%rootFunction(densityContrastsRoot         )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
    end if
    ! Make the self, node, and basic component available to the root finding routine.
    densityContrastsSelf  => self
    densityContrastsNode  => node
    densityContrastsBasic => node%basic()
    do i=1,self%countDensityContrasts
       densityContrastsMeanTarget=self  %densityContrasts(          i                                           )
       radius                    =finder%find            (rootGuess=self%darkMatterHaloScale_%virialRadius(node))
       enclosedMass              =Galactic_Structure_Enclosed_Mass(                                     &
            &                                                                         node            , &
            &                                                                         radius          , &
            &                                                      componentType=     componentTypeAll, &
            &                                                      massType     =self%massTypeSelected  &
            &                                                     )
       densityContrastsExtract((i-1)*2+1:(i-1)*2+2)=[radius,enclosedMass]
    end do
    return
  end function densityContrastsExtract

  function densityContrastsNames(self,time)
    !% Return the names of the {\normalfont \ttfamily densityContrasts} properties.
    implicit none
    type            (varying_string                       ), dimension(:) , allocatable :: densityContrastsNames
    class           (nodePropertyExtractorDensityContrasts), intent(inout)              :: self
    double precision                                       , intent(in   )              :: time
    integer                                                                             :: i
    character       (len=32                               )                             :: name                 , label
    !$GLC attributes unused :: time

    allocate(densityContrastsNames(self%elementCount_))
    do i=1,self%countDensityContrasts
       write (label,'(f9.4)') self%densityContrasts(i)
       write (name ,'(a,a)' ) 'nodeRadius',trim(adjustl(label))
       densityContrastsNames((i-1)*2+1)=trim(adjustl(name))
       write (name ,'(a,a)' ) 'nodeMass'  ,trim(adjustl(label))
       densityContrastsNames((i-1)*2+2)=trim(adjustl(name))
    end do
    return
  end function densityContrastsNames

  function densityContrastsDescriptions(self,time)
    !% Return descriptions of the {\normalfont \ttfamily densityContrasts} property.
    implicit none
    type            (varying_string                       ), dimension(:) , allocatable :: densityContrastsDescriptions
    class           (nodePropertyExtractorDensityContrasts), intent(inout)              :: self
    double precision                                       , intent(in   )              :: time
    integer                                                                             :: i
    character       (len=64                               )                             :: description
    !$GLC attributes unused :: time

    allocate(densityContrastsDescriptions(self%elementCount_))
    do i=1,self%countDensityContrasts
       write (description,'(a,f5.1,a)') 'Radius enclosing a density contrast of ',self%densityContrasts(i),' [Mpc].'
       densityContrastsDescriptions((i-1)*2+1)=trim(adjustl(description))
       write (description,'(a,f5.1,a)') 'Mass within a density contrast of '     ,self%densityContrasts(i),' [M☉].'
       densityContrastsDescriptions((i-1)*2+2)=trim(adjustl(description))
    end do
    return
  end function densityContrastsDescriptions

  function densityContrastsUnitsInSI(self,time)
    !% Return the units of the {\normalfont \ttfamily densityContrasts} properties in the SI system.
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                       , allocatable  , dimension(:) :: densityContrastsUnitsInSI
    class           (nodePropertyExtractorDensityContrasts), intent(inout)               :: self
    double precision                                       , intent(in   )               :: time
    integer                                                                              :: i
    !$GLC attributes unused :: time

    allocate(densityContrastsUnitsInSI(self%elementCount_))
    do i=1,self%countDensityContrasts
       densityContrastsUnitsInSI((i-1)*2+1:(i-1)*2+2)=[megaParsec,massSolar]
    end do
    return
  end function densityContrastsUnitsInSI

  integer function densityContrastsType(self)
    !% Return the type of the {\normalfont \ttfamily densityContrasts} properties.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorDensityContrasts), intent(inout) :: self
    !$GLC attributes unused :: self

    densityContrastsType=outputAnalysisPropertyTypeLinear
    return
  end function densityContrastsType

  double precision function densityContrastsRoot(radius)
    !% Root function used in finding the radius that encloses a given density contrast.
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options        , only : componentTypeAll
    use :: Numerical_Constants_Math          , only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: enclosedMass

    enclosedMass        =Galactic_Structure_Enclosed_Mass(                                                         &
         &                                                                                   densityContrastsNode, &
         &                                                                                   radius              , &
         &                                                componentType=                     componentTypeAll    , &
         &                                                massType     =densityContrastsSelf%massTypeSelected      &
         &                                               )
    densityContrastsRoot=+3.0d0                                                                                       &
         &               *enclosedMass                                                                                &
         &               /4.0d0                                                                                       &
         &               /Pi                                                                                          &
         &               /radius**3                                                                                   &
         &               /(                                                                                           &
         &                 +densityContrastsSelf%densityReference                                                     &
         &                 /densityContrastsSelf%cosmologyFunctions_%expansionFactor(densityContrastsBasic%time())**3 &
         &               )                                                                                            &
         &               -densityContrastsMeanTarget
    return
  end function densityContrastsRoot
