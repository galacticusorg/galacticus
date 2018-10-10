!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implements a merger mass movements class using the \cite{baugh_can_2005} model.
  
  !# <mergerMassMovements name="mergerMassMovementsBaugh2005">
  !#  <description>A merger mass movements class which uses a simple calculation.</description>
  !# </mergerMassMovements>
  type, extends(mergerMassMovementsClass) :: mergerMassMovementsBaugh2005
     !% A merger mass movements class which uses the \cite{baugh_can_2005} calculation.
     private
     double precision :: massRatioMajorMerger     , ratioMassBurst, &
          &              fractionGasCriticalBurst
     integer          :: destinationGasMinorMerger
   contains
     procedure :: get => baugh2005Get
  end type mergerMassMovementsBaugh2005

  interface mergerMassMovementsBaugh2005
     !% Constructors for the {\normalfont \ttfamily baugh2005} merger mass movements class.
     module procedure baugh2005ConstructorParameters
     module procedure baugh2005ConstructorInternal
  end interface mergerMassMovementsBaugh2005

contains

  function baugh2005ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily baugh2005} merger mass movements class which takes a parameter list as input.
    use Input_Parameters
    implicit none
    type            (mergerMassMovementsBaugh2005)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: massRatioMajorMerger     , ratioMassBurst, &
         &                                                           fractionGasCriticalBurst
    type            (varying_string              )                :: destinationGasMinorMerger

    !# <inputParameter>
    !#   <name>massRatioMajorMerger</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.25d0</defaultValue>
    !#   <description>The mass ratio above which mergers are considered to be ``major''.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>ratioMassBurst</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.05d0</defaultValue>
    !#   <description>The mass ratio above which mergers are considered to trigger a burst.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>fractionGasCriticalBurst</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.75d0</defaultValue>
    !#   <description>The host gas fraction above which mergers are considered to trigger a burst.</description>
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
    self=mergerMassMovementsBaugh2005(massRatioMajorMerger,enumerationDestinationMergerEncode(char(destinationGasMinorMerger),includesPrefix=.false.),ratioMassBurst,fractionGasCriticalBurst)
    !# <inputParametersValidate source="parameters"/>
    return
  end function baugh2005ConstructorParameters

  function baugh2005ConstructorInternal(massRatioMajorMerger,destinationGasMinorMerger,ratioMassBurst,fractionGasCriticalBurst) result(self)
    !% Internal constructor for the {\normalfont \ttfamily baugh2005} merger mass movements.
    implicit none
    type            (mergerMassMovementsBaugh2005)                :: self
    double precision                              , intent(in   ) :: massRatioMajorMerger     , ratioMassBurst, &
         &                                                           fractionGasCriticalBurst
    integer                                       , intent(in   ) :: destinationGasMinorMerger
    !# <constructorAssign variables="massRatioMajorMerger, destinationGasMinorMerger, ratioMassBurst, fractionGasCriticalBurst"/>
    
    return
  end function baugh2005ConstructorInternal

  subroutine baugh2005Get(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !% Determine how different mass components should be redistributed as the result of a merger according to the model of
    !% \cite{baugh_can_2005}.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    class           (mergerMassMovementsBaugh2005), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    integer                                       , intent(  out) :: destinationGasSatellite, destinationGasHost       , &
         &                                                           destinationStarsHost   , destinationStarsSatellite
    logical                                       , intent(  out) :: mergerIsMajor
    type            (treeNode                    ), pointer       :: nodeHost
    double precision                                              :: massHost               , massSatellite            , &
         &                                                           massSpheroidHost       , massGasHost
    logical                                                       :: triggersBurst
    
    nodeHost         => node%mergesWith()
    massSatellite    =  Galactic_Structure_Enclosed_Mass(node                                        ,massType=massTypeGalactic)
    massHost         =  Galactic_Structure_Enclosed_Mass(nodeHost                                    ,massType=massTypeGalactic)
    massGasHost      =  Galactic_Structure_Enclosed_Mass(nodeHost                                    ,massType=massTypeGaseous )
    massSpheroidHost =  Galactic_Structure_Enclosed_Mass(nodeHost,componentType=componentTypeSpheroid,massType=massTypeGalactic)
    mergerIsMajor    =  massSatellite >= self%massRatioMajorMerger*massHost

    triggersBurst=mergerIsMajor                                               &
         &         .or.                                                       &
         &        (                                                           &
         &         massSpheroidHost <  self%ratioMassBurst          *massHost &
         &          .and.                                                     &
         &         massGasHost      >= self%fractionGasCriticalBurst*massHost &
         &        )
    if (mergerIsMajor) then
       destinationGasSatellite     =    destinationMergerSpheroid
       destinationStarsSatellite   =    destinationMergerSpheroid
       destinationGasHost          =    destinationMergerSpheroid
       destinationStarsHost        =    destinationMergerSpheroid
    else
       if (triggersBurst) then
          destinationGasSatellite  =    destinationMergerSpheroid
          destinationStarsSatellite=    destinationMergerSpheroid
          destinationGasHost       =    destinationMergerSpheroid
          destinationStarsHost     =    destinationMergerUnmoved
       else
          destinationGasSatellite  =self%destinationGasMinorMerger
          destinationStarsSatellite=    destinationMergerSpheroid
          destinationGasHost       =    destinationMergerUnmoved
          destinationStarsHost     =    destinationMergerUnmoved
       end if
    end if
    return
  end subroutine baugh2005Get
