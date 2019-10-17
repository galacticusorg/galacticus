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

  !# <mergerMassMovements name="mergerMassMovementsSimple">
  !#  <description>A merger mass movements class which uses a simple calculation.</description>
  !# </mergerMassMovements>
  type, extends(mergerMassMovementsClass) :: mergerMassMovementsSimple
     !% A merger mass movements class which uses a simple calculation.
     private
     double precision :: massRatioMajorMerger
     integer          :: destinationGasMinorMerger
   contains
     procedure :: get => simpleGet
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

    return
  end function simpleConstructorInternal

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

    nodeHost      => node%mergesWith()
    massSatellite =  Galactic_Structure_Enclosed_Mass(node    ,massType=massTypeGalactic)
    massHost      =  Galactic_Structure_Enclosed_Mass(nodeHost,massType=massTypeGalactic)
    mergerIsMajor =  massSatellite >= self%massRatioMajorMerger*massHost
    if (mergerIsMajor) then
       destinationGasSatellite  =    destinationMergerSpheroid
       destinationStarsSatellite=    destinationMergerSpheroid
       destinationGasHost       =    destinationMergerSpheroid
       destinationStarsHost     =    destinationMergerSpheroid
    else
       destinationGasSatellite  =self%destinationGasMinorMerger
       destinationStarsSatellite=    destinationMergerSpheroid
       destinationGasHost       =    destinationMergerUnmoved
       destinationStarsHost     =    destinationMergerUnmoved
    end if
    return
  end subroutine simpleGet
