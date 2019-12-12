!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a merger mass movements class which uses a simple calculation.

  use :: Kind_Numbers, only : kind_int8

  !# <mergerMassMovements name="mergerMassMovementsSimple">
  !#  <description>A merger mass movements class which uses a simple calculation.</description>
  !# </mergerMassMovements>
  type, extends(mergerMassMovementsClass) :: mergerMassMovementsSimple
     !% A merger mass movements class which uses a simple calculation.
     private
     double precision                 :: massRatioMajorMerger
     integer                          :: destinationGasMinorMerger
     integer         (kind=kind_int8) :: lastUniqueID
     integer                          :: destinationGasSatellite  , destinationStarsSatellite, &
          &                              destinationGasHost       , destinationStarsHost
     logical                          :: mergerIsMajor            , movementsCalculated
   contains
     final     ::             simpleDestructor
     procedure :: autoHook => simpleAutoHook
     procedure :: get      => simpleGet
  end type mergerMassMovementsSimple

  interface mergerMassMovementsSimple
     !% Constructors for the {\normalfont \ttfamily simple} merger mass movements class.
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface mergerMassMovementsSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily simple} merger mass movements class which takes a parameter list as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerMassMovementsSimple)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: massRatioMajorMerger
    type            (varying_string           )                :: destinationGasMinorMerger

    !# <inputParameter>
    !#   <name>massRatioMajorMerger</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.25d0</defaultValue>
    !#   <description>The mass ratio above which mergers are considered to be ``major''.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>destinationGasMinorMerger</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>var_str('spheroid')</defaultValue>
    !#   <description>The component to which satellite galaxy gas moves to as a result of a minor merger.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    self=mergerMassMovementsSimple(massRatioMajorMerger,enumerationDestinationMergerEncode(char(destinationGasMinorMerger),includesPrefix=.false.))
    !# <inputParametersValidate source="parameters"/>
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(massRatioMajorMerger,destinationGasMinorMerger) result(self)
    !% Internal constructor for the {\normalfont \ttfamily simple} merger mass movements class.
    implicit none
    type            (mergerMassMovementsSimple)                :: self
    double precision                           , intent(in   ) :: massRatioMajorMerger
    integer                                    , intent(in   ) :: destinationGasMinorMerger
    !# <constructorAssign variables="massRatioMajorMerger, destinationGasMinorMerger"/>

    self%lastUniqueID             =-huge(0_kind_int8)
    self%destinationGasSatellite  =-huge(0          )
    self%destinationStarsSatellite=-huge(0          )
    self%destinationGasHost       =-huge(0          )
    self%destinationStarsHost     =-huge(0          )
    self%mergerIsMajor            =.false.
    self%movementsCalculated      =.false.
    return
  end function simpleConstructorInternal

  subroutine simpleAutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, satelliteMergerEvent, openMPThreadBindingAllLevels
    implicit none
    class(mergerMassMovementsSimple), intent(inout) :: self

    call calculationResetEvent%attach(self,simpleCalculationReset,openMPThreadBindingAllLevels                                             )
    call satelliteMergerEvent %attach(self,simpleGetHook         ,openMPThreadBindingAllLevels,label='remnantStructure:massMovementsSimple')
    return
  end subroutine simpleAutoHook

  subroutine simpleDestructor(self)
    !% Destructor for the {\normalfont \ttfamily simple} satellite merger mass movements class
    use :: Events_Hooks, only : calculationResetEvent, satelliteMergerEvent
    implicit none
    type(mergerMassMovementsSimple), intent(inout) :: self

    call calculationResetEvent%detach(self,simpleCalculationReset)
    call satelliteMergerEvent %detach(self,simpleGetHook         )
    return
  end subroutine simpleDestructor

  subroutine simpleCalculationReset(self,node)
    !% Reset the dark matter profile calculation.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(*       ), intent(inout) :: self
    type (treeNode), intent(inout) :: node

    select type (self)
    class is (mergerMassMovementsSimple)
       self%movementsCalculated=.false.
       self%lastUniqueID       =node%uniqueID()
       class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine simpleCalculationReset

  subroutine simpleGetHook(self,node)
    !% Hookable wrapper around the get function.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class  (*       ), intent(inout)         :: self
    type   (treeNode), intent(inout), target :: node
    integer                                  :: destinationGasSatellite, destinationGasHost       , &
         &                                      destinationStarsHost   , destinationStarsSatellite
    logical                                  :: mergerIsMajor

    select type (self)
    type is (mergerMassMovementsSimple)
       call self%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine simpleGetHook

  subroutine simpleGet(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !% Determine where stars and gas move as the result of a merger event using a simple algorithm.
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options        , only : massTypeGalactic
    implicit none
    class           (mergerMassMovementsSimple), intent(inout) :: self
    type            (treeNode                 ), intent(inout) :: node
    integer                                    , intent(  out) :: destinationGasSatellite, destinationGasHost       , &
         &                                                        destinationStarsHost   , destinationStarsSatellite
    logical                                    , intent(  out) :: mergerIsMajor
    type            (treeNode                 ), pointer       :: nodeHost
    double precision                                           :: massHost               , massSatellite

    ! The calculation of how mass moves as a result of the merger is computed when first needed and then stored. This ensures that
    ! the results are determined by the properties of the merge target prior to any modification that will occur as node
    ! components are modified in response to the merger.
    if (node%uniqueID() /= self%lastUniqueID) call simpleCalculationReset(self,node)
    if (.not.self%movementsCalculated) then
       self%movementsCalculated =  .true.
       nodeHost                 => node%mergesWith()
       massSatellite            =  Galactic_Structure_Enclosed_Mass(node    ,massType=massTypeGalactic)
       massHost                 =  Galactic_Structure_Enclosed_Mass(nodeHost,massType=massTypeGalactic)
       self%mergerIsMajor       =  massSatellite >= self%massRatioMajorMerger*massHost
       if (self%mergerIsMajor) then
          self%destinationGasSatellite  =    destinationMergerSpheroid
          self%destinationStarsSatellite=    destinationMergerSpheroid
          self%destinationGasHost       =    destinationMergerSpheroid
          self%destinationStarsHost     =    destinationMergerSpheroid
       else
          self%destinationGasSatellite  =self%destinationGasMinorMerger
          self%destinationStarsSatellite=    destinationMergerSpheroid
          self%destinationGasHost       =    destinationMergerUnmoved
          self%destinationStarsHost     =    destinationMergerUnmoved
       end if
    end if
    mergerIsMajor            =self%mergerIsMajor
    destinationGasSatellite  =self%destinationGasSatellite
    destinationStarsSatellite=self%destinationStarsSatellite
    destinationGasHost       =self%destinationGasHost
    destinationStarsHost     =self%destinationStarsHost
    return
  end subroutine simpleGet
