!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a merger mass movements class which uses a simple calculation.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <mergerMassMovements name="mergerMassMovementsSimple">
   <description>
    A merger mass movements class which implements mass movements according to:
    \begin{itemize}
     \item If $M_\mathrm{satellite} &gt; f_\mathrm{major} M_\mathrm{central}$ then all mass from both satellite and central
     galaxies moves to the spheroid \gls{component} of the central galaxy;
     \item Otherwise: Gas from the satellite moves to the \gls{component} of the central specified by the {\normalfont
     \ttfamily [minorMergerGasMovesTo]} parameter (either ``{\normalfont \ttfamily disk}'' or ``{\normalfont \ttfamily
     spheroid}''), stars from the satellite moves to the spheroid of the central and mass in the central does not move.
    \end{itemize}
    Here, $f_\mathrm{major}=${\normalfont \ttfamily [majorMergerMassRatio]} is the mass ratio above which a merger is
    considered to be ``major''.
   </description>
  </mergerMassMovements>
  !!]
  type, extends(mergerMassMovementsClass) :: mergerMassMovementsSimple
     !!{
     A merger mass movements class which uses a simple calculation.
     !!}
     private
     double precision                                   :: massRatioMajorMerger
     type            (enumerationDestinationMergerType) :: destinationGasMinorMerger, destinationStarsMinorMerger
     integer         (kind=kind_int8                  ) :: lastUniqueID
     type            (enumerationDestinationMergerType) :: destinationGasSatellite  , destinationStarsSatellite  , &
          &                                                destinationGasHost       , destinationStarsHost
     logical                                            :: mergerIsMajor            , movementsCalculated
   contains
     final     ::             simpleDestructor
     procedure :: autoHook => simpleAutoHook
     procedure :: get      => simpleGet
  end type mergerMassMovementsSimple

  interface mergerMassMovementsSimple
     !!{
     Constructors for the \refClass{mergerMassMovementsSimple} merger mass movements class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface mergerMassMovementsSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerMassMovementsSimple} merger mass movements class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerMassMovementsSimple)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: massRatioMajorMerger
    type            (varying_string           )                :: destinationGasMinorMerger, destinationStarsMinorMerger

    !![
    <inputParameter>
      <name>massRatioMajorMerger</name>
      <defaultValue>0.25d0</defaultValue>
      <description>The mass ratio above which mergers are considered to be ``major''.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>destinationGasMinorMerger</name>
      <defaultValue>var_str('spheroid')</defaultValue>
      <description>The component to which satellite galaxy gas moves to as a result of a minor merger.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>destinationStarsMinorMerger</name>
      <defaultValue>var_str('spheroid')</defaultValue>
      <description>The component to which satellite galaxy stars move to as a result of a minor merger.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerMassMovementsSimple(massRatioMajorMerger,enumerationDestinationMergerEncode(char(destinationGasMinorMerger),includesPrefix=.false.),enumerationDestinationMergerEncode(char(destinationStarsMinorMerger),includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(massRatioMajorMerger,destinationGasMinorMerger,destinationStarsMinorMerger) result(self)
    !!{
    Internal constructor for the \refClass{mergerMassMovementsSimple} merger mass movements class.
    !!}
    implicit none
    type            (mergerMassMovementsSimple       )                        :: self
    double precision                                  , intent(in   )         :: massRatioMajorMerger
    type            (enumerationDestinationMergerType), intent(in   )         :: destinationGasMinorMerger, destinationStarsMinorMerger
    !![
    <constructorAssign variables="massRatioMajorMerger, destinationGasMinorMerger, destinationStarsMinorMerger"/>
    !!]

    self%lastUniqueID             =-huge(0_kind_int8)
    self%destinationGasSatellite  =destinationMergerUnmoved
    self%destinationStarsSatellite=destinationMergerUnmoved
    self%destinationGasHost       =destinationMergerUnmoved
    self%destinationStarsHost     =destinationMergerUnmoved
    self%mergerIsMajor            =.false.
    self%movementsCalculated      =.false.
    return
  end function simpleConstructorInternal

  subroutine simpleAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels, satelliteMergerEvent
    implicit none
    class(mergerMassMovementsSimple), intent(inout) :: self

    call calculationResetEvent%attach(self,simpleCalculationReset,openMPThreadBindingAllLevels,label='remnantStructure:massMovementsSimple')
    call satelliteMergerEvent %attach(self,simpleGetHook         ,openMPThreadBindingAllLevels,label='remnantStructure:massMovementsSimple')
    return
  end subroutine simpleAutoHook

  subroutine simpleDestructor(self)
    !!{
    Destructor for the \refClass{mergerMassMovementsSimple} satellite merger mass movements class
    !!}
    use :: Events_Hooks, only : calculationResetEvent, satelliteMergerEvent
    implicit none
    type(mergerMassMovementsSimple), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,simpleCalculationReset)) call calculationResetEvent%detach(self,simpleCalculationReset)
    if (satelliteMergerEvent %isAttached(self,simpleGetHook         )) call satelliteMergerEvent %detach(self,simpleGetHook         )
    return
  end subroutine simpleDestructor

  subroutine simpleCalculationReset(self,node,uniqueID)
    !!{
    Reset the dark matter profile calculation.
    !!}
    use :: Error       , only : Error_Report
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (*        ), intent(inout) :: self
    type   (treeNode ), intent(inout) :: node
    integer(kind_int8), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    select type (self)
    class is (mergerMassMovementsSimple)
       self%movementsCalculated=.false.
       self%lastUniqueID       =uniqueID
       class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine simpleCalculationReset

  subroutine simpleGetHook(self,node)
    !!{
    Hookable wrapper around the get function.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (*                               ), intent(inout)         :: self
    type   (treeNode                        ), intent(inout), target :: node
    type   (enumerationDestinationMergerType)                        :: destinationGasSatellite, destinationGasHost       , &
         &                                                              destinationStarsHost   , destinationStarsSatellite
    logical                                                          :: mergerIsMajor

    select type (self)
    type is (mergerMassMovementsSimple)
       call self%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine simpleGetHook

  subroutine simpleGet(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !!{
    Determine where stars and gas move as the result of a merger event using a simple algorithm.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk    , componentTypeSpheroid, massTypeGalactic
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (mergerMassMovementsSimple       ), intent(inout)         :: self
    type            (treeNode                        ), intent(inout), target :: node
    type            (enumerationDestinationMergerType), intent(  out)         :: destinationGasSatellite  , destinationGasHost       , &
         &                                                                       destinationStarsHost     , destinationStarsSatellite
    logical                                           , intent(  out)         :: mergerIsMajor
    type            (treeNode                        ), pointer               :: nodeHost                 , nodeMajor
    class           (massDistributionClass           ), pointer               :: massDistributionSatellite, massDistributionHost     , &
         &                                                                       massDistributionDisk     , massDistributionSpheroid
    double precision                                                          :: massHost                 , massSatellite            , &
         &                                                                       massSpheroid             , massDisk
    type            (enumerationDestinationMergerType)                        :: destinationDominant

    ! The calculation of how mass moves as a result of the merger is computed when first needed and then stored. This ensures that
    ! the results are determined by the properties of the merge target prior to any modification that will occur as node
    ! components are modified in response to the merger.
    if (node%uniqueID() /= self%lastUniqueID) call simpleCalculationReset(self,node,node%uniqueID())
    if (.not.self%movementsCalculated) then
       self%movementsCalculated  =  .true.
       nodeHost                  => node                     %mergesWith      (                         )
       massDistributionHost      => nodeHost                 %massDistribution(massType=massTypeGalactic)
       massDistributionSatellite => node                     %massDistribution(massType=massTypeGalactic)
       massSatellite             =  massDistributionSatellite%massTotal       (                         )
       massHost                  =  massDistributionHost     %massTotal       (                         )
       self%mergerIsMajor        =  massSatellite > 0.0d0 .and. massHost > 0.0d0 .and. min(massSatellite,massHost) >= self%massRatioMajorMerger*max(massSatellite,massHost)
       !![
       <objectDestructor name="massDistributionHost"     />
       <objectDestructor name="massDistributionSatellite"/>
       !!]
       if (self%mergerIsMajor) then
          self%destinationGasSatellite  =     destinationMergerSpheroid
          self%destinationStarsSatellite=     destinationMergerSpheroid
          self%destinationGasHost       =     destinationMergerSpheroid
          self%destinationStarsHost     =     destinationMergerSpheroid
       else
          destinationDominant=destinationMergerUnmoved
          if (self%destinationGasMinorMerger == destinationMergerDominant .or. self%destinationStarsMinorMerger == destinationMergerDominant) then
             if (massSatellite < massHost) then
                nodeMajor => nodeHost
             else
                nodeMajor => node
             end if
             massDistributionDisk     => nodeMajor               %massDistribution(massType=massTypeGalactic,componentType=componentTypeDisk    )
             massDistributionSpheroid => nodeMajor               %massDistribution(massType=massTypeGalactic,componentType=componentTypeSpheroid)
             massDisk                 =  massDistributionDisk    %massTotal       (                                                             )
             massSpheroid             =  massDistributionSpheroid%massTotal       (                                                             )
             !![
	     <objectDestructor name="massDistributionDisk"    />
	     <objectDestructor name="massDistributionSpheroid"/>
	     !!]
             if (massDisk > massSpheroid) then
                destinationDominant=destinationMergerDisk
             else
                destinationDominant=destinationMergerSpheroid
             end if
          end if
          if (self%destinationGasMinorMerger   == destinationMergerDominant) then
             self%destinationGasSatellite  =     destinationDominant
          else
             self%destinationGasSatellite  =self%destinationGasMinorMerger
          end if
          if (self%destinationStarsMinorMerger == destinationMergerDominant) then
             self%destinationStarsSatellite=     destinationDominant
          else
             self%destinationStarsSatellite=self%destinationStarsMinorMerger
          end if
          self%destinationGasHost       =     destinationMergerUnmoved
          self%destinationStarsHost     =     destinationMergerUnmoved
       end if
    end if
    mergerIsMajor            =self%mergerIsMajor
    destinationGasSatellite  =self%destinationGasSatellite
    destinationStarsSatellite=self%destinationStarsSatellite
    destinationGasHost       =self%destinationGasHost
    destinationStarsHost     =self%destinationStarsHost
    return
  end subroutine simpleGet
