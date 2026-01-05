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
  Implements a merger mass movements class using the \cite{baugh_can_2005} model.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <mergerMassMovements name="mergerMassMovementsBaugh2005">
   <description>
    A merger mass movements class which implements mass movements according to:
    \begin{itemize}
     \item If $M_\mathrm{satellite} &gt; f_\mathrm{major} M_\mathrm{central}$ then all mass from both satellite and central
     galaxies moves to the spheroid \gls{component} of the central galaxy;
     \item Otherwise:
     \begin{itemize}
      \item If $M_\mathrm{central, spheroid} &lt; f_\mathrm{burst} M_\mathrm{central}$ and the gas fraction in the host equals or
      exceeds $f_\mathrm{gas,crit}$ then all gas is moved to the host spheroid, while the host stellar disk remains in place.
      \item Otherwise, gas from the satellite moves to the \gls{component} of the central specified by the {\normalfont
      \ttfamily [destinationGasMinorMerger]} parameter (either ``{\normalfont \ttfamily disk}'' or ``{\normalfont \ttfamily
      spheroid}''), stars from the satellite moves to the spheroid of the central and mass in the central does not move.
     \end{itemize}
    \end{itemize}
    Here, $f_\mathrm{major}=${\normalfont \ttfamily [massRatioMajorMerger]} is the mass ratio above which a merger is
    considered to be ``major'', while $f_\mathrm{burst}=${\normalfont \ttfamily [ratioMassBurst]} and
    $f_\mathrm{gas,crit}=${\normalfont \ttfamily [fractionGasCriticalBurst]}.
   </description>
  </mergerMassMovements>
  !!]
  type, extends(mergerMassMovementsClass) :: mergerMassMovementsBaugh2005
     !!{
     A merger mass movements class which uses the \cite{baugh_can_2005} calculation.
     !!}
     private
     double precision                                   :: massRatioMajorMerger     , ratioMassBurst           , &
          &                                                fractionGasCriticalBurst
     type            (enumerationDestinationMergerType) :: destinationGasMinorMerger
     integer         (kind=kind_int8                  ) :: lastUniqueID
     type            (enumerationDestinationMergerType) :: destinationGasSatellite  , destinationStarsSatellite, &
          &                                                destinationGasHost       , destinationStarsHost
     logical                                            :: mergerIsMajor            , movementsCalculated
   contains
     final     ::             baugh2005Destructor
     procedure :: autoHook => baugh2005AutoHook
     procedure :: get      => baugh2005Get
  end type mergerMassMovementsBaugh2005

  interface mergerMassMovementsBaugh2005
     !!{
     Constructors for the \refClass{mergerMassMovementsBaugh2005} merger mass movements class.
     !!}
     module procedure baugh2005ConstructorParameters
     module procedure baugh2005ConstructorInternal
  end interface mergerMassMovementsBaugh2005

contains

  function baugh2005ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerMassMovementsBaugh2005} merger mass movements class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerMassMovementsBaugh2005)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: massRatioMajorMerger     , ratioMassBurst, &
         &                                                           fractionGasCriticalBurst
    type            (varying_string              )                :: destinationGasMinorMerger

    !![
    <inputParameter>
      <name>massRatioMajorMerger</name>
      <defaultValue>0.25d0</defaultValue>
      <description>The mass ratio above which mergers are considered to be ``major''.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>ratioMassBurst</name>
      <defaultValue>0.05d0</defaultValue>
      <description>The mass ratio above which mergers are considered to trigger a burst.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fractionGasCriticalBurst</name>
      <defaultValue>0.75d0</defaultValue>
      <description>The host gas fraction above which mergers are considered to trigger a burst.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>destinationGasMinorMerger</name>
      <defaultValue>var_str('spheroid')</defaultValue>
      <description>The component to which satellite galaxy gas moves to as a result of a minor merger.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerMassMovementsBaugh2005(massRatioMajorMerger,enumerationDestinationMergerEncode(char(destinationGasMinorMerger),includesPrefix=.false.),ratioMassBurst,fractionGasCriticalBurst)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function baugh2005ConstructorParameters

  function baugh2005ConstructorInternal(massRatioMajorMerger,destinationGasMinorMerger,ratioMassBurst,fractionGasCriticalBurst) result(self)
    !!{
    Internal constructor for the \refClass{mergerMassMovementsBaugh2005} merger mass movements.
    !!}
    implicit none
    type            (mergerMassMovementsBaugh2005    )                :: self
    double precision                                  , intent(in   ) :: massRatioMajorMerger     , ratioMassBurst, &
         &                                                               fractionGasCriticalBurst
    type            (enumerationDestinationMergerType), intent(in   ) :: destinationGasMinorMerger
    !![
    <constructorAssign variables="massRatioMajorMerger, destinationGasMinorMerger, ratioMassBurst, fractionGasCriticalBurst"/>
    !!]

    self%lastUniqueID             =-huge(0_kind_int8)
    self%destinationGasSatellite  =destinationMergerUnmoved
    self%destinationStarsSatellite=destinationMergerUnmoved
    self%destinationGasHost       =destinationMergerUnmoved
    self%destinationStarsHost     =destinationMergerUnmoved
    self%mergerIsMajor            =.false.
    self%movementsCalculated      =.false.
    return
  end function baugh2005ConstructorInternal

  subroutine baugh2005AutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels, satelliteMergerEvent
    implicit none
    class(mergerMassMovementsBaugh2005), intent(inout) :: self

    call calculationResetEvent%attach(self,baugh2005CalculationReset,openMPThreadBindingAllLevels,label='remnantStructure:massMovementsBaugh2005')
    call satelliteMergerEvent %attach(self,baugh2005GetHook         ,openMPThreadBindingAllLevels,label='remnantStructure:massMovementsBaugh2005')
    return
  end subroutine baugh2005AutoHook

  subroutine baugh2005Destructor(self)
    !!{
    Destructor for the \refClass{mergerMassMovementsBaugh2005} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, satelliteMergerEvent
    implicit none
    type(mergerMassMovementsBaugh2005), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,baugh2005CalculationReset)) call calculationResetEvent%detach(self,baugh2005CalculationReset)
    if (satelliteMergerEvent %isAttached(self,baugh2005GetHook         )) call satelliteMergerEvent %detach(self,baugh2005GetHook         )
    return
  end subroutine baugh2005Destructor

  subroutine baugh2005CalculationReset(self,node,uniqueID)
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
    class is (mergerMassMovementsBaugh2005)
       self%movementsCalculated=.false.
       self%lastUniqueID       =uniqueID
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine baugh2005CalculationReset

  subroutine baugh2005GetHook(self,node)
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
    type is (mergerMassMovementsBaugh2005)
       call self%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine baugh2005GetHook

  subroutine baugh2005Get(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !!{
    Determine how different mass components should be redistributed as the result of a merger according to the model of
    \cite{baugh_can_2005}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeSpheroid, componentTypeDisk, massTypeGalactic, massTypeGaseous
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (mergerMassMovementsBaugh2005    ), intent(inout)         :: self
    type            (treeNode                        ), intent(inout), target :: node
    type            (enumerationDestinationMergerType), intent(  out)         :: destinationGasSatellite     , destinationGasHost             , &
         &                                                                       destinationStarsHost        , destinationStarsSatellite
    logical                                           , intent(  out)         :: mergerIsMajor
    type            (treeNode                        ), pointer               :: nodeHost
    class           (massDistributionClass           ), pointer               :: massDistributionSatellite   , massDistributionHost           , &
         &                                                                       massDistributionHostDiskGas , massDistributionHostSpheroidGas, &
         &                                                                       massDistributionHostSpheroid
    double precision                                                          :: massHost                    , massSatellite                  , &
         &                                                                       massSpheroidHost            , massGasHost
    logical                                                                   :: triggersBurst
    
    ! The calculation of how mass moves as a result of the merger is computed when first needed and then stored. This ensures that
    ! the results are determined by the properties of the merge target prior to any modification that will occur as node
    ! components are modified in response to the merger.
    if (node%uniqueID() /= self%lastUniqueID) call baugh2005CalculationReset(self,node,node%uniqueID())
    if (.not.self%movementsCalculated) then
       self%movementsCalculated =  .true.
       nodeHost                        =>  node                           %mergesWith      (                                                             )
       massDistributionSatellite       =>  node                           %massDistribution(                                    massType=massTypeGalactic) 
       massDistributionHost            =>  nodeHost                       %massDistribution(                                    massType=massTypeGalactic)
       massDistributionHostSpheroid    =>  nodeHost                       %massDistribution(componentType=componentTypeSpheroid,massType=massTypeGalactic)
       massDistributionHostDiskGas     =>  nodeHost                       %massDistribution(componentType=componentTypeDisk    ,massType=massTypeGaseous )
       massDistributionHostSpheroidGas =>  nodeHost                       %massDistribution(componentType=componentTypeSpheroid,massType=massTypeGaseous )
       massSatellite                   =  +massDistributionSatellite      %massTotal       (                                                             )
       massHost                        =  +massDistributionHost           %massTotal       (                                                             )
       massGasHost                     =  +massDistributionHostDiskGas    %massTotal       (                                                             ) &
            &                             +massDistributionHostSpheroidGas%massTotal       (                                                             )
       massSpheroidHost                =  +massDistributionHostSpheroid   %massTotal       (                                                             )
       self%mergerIsMajor              =    massSatellite    >= self%massRatioMajorMerger    *massHost
       triggersBurst                   =                        self%mergerIsMajor                     &
            &                             .or.                                                         &
            &                              (                                                           &
            &                               massSpheroidHost <  self%ratioMassBurst          *massHost &
            &                                .and.                                                     &
            &                               massGasHost      >= self%fractionGasCriticalBurst*massHost &
            &                              )
       !![
       <objectDestructor name="massDistributionSatellite"      />
       <objectDestructor name="massDistributionHost"           />
       <objectDestructor name="massDistributionHostDiskGas"    />
       <objectDestructor name="massDistributionHostSpheroidGas"/>
       <objectDestructor name="massDistributionHostSpheroid"   />
       !!]
        if (self%mergerIsMajor) then
          self%destinationGasSatellite     =    destinationMergerSpheroid
          self%destinationStarsSatellite   =    destinationMergerSpheroid
          self%destinationGasHost          =    destinationMergerSpheroid
          self%destinationStarsHost        =    destinationMergerSpheroid
       else
          if (triggersBurst) then
             self%destinationGasSatellite  =    destinationMergerSpheroid
             self%destinationStarsSatellite=    destinationMergerSpheroid
             self%destinationGasHost       =    destinationMergerSpheroid
             self%destinationStarsHost     =    destinationMergerUnmoved
          else
             self%destinationGasSatellite  =self%destinationGasMinorMerger
             self%destinationStarsSatellite=    destinationMergerSpheroid
             self%destinationGasHost       =    destinationMergerUnmoved
             self%destinationStarsHost     =    destinationMergerUnmoved
          end if
       end if
    end if
    mergerIsMajor            =self%mergerIsMajor
    destinationGasSatellite  =self%destinationGasSatellite
    destinationStarsSatellite=self%destinationStarsSatellite
    destinationGasHost       =self%destinationGasHost
    destinationStarsHost     =self%destinationStarsHost
    return
  end subroutine baugh2005Get
