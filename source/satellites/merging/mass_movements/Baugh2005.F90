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

  !+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a merger mass movements class using the :cite:t:`baugh_can_2005` model.
  !!}

  ! Enumeration for the criterion used to decide whether a minor merger triggers a burst.
  !![
  <enumeration docformat="rst">
   <name>baugh2005BurstCriterion</name>
   <description>
   Specifies the criterion used by the :galacticus-class:`mergerMassMovementsBaugh2005` class to decide whether a minor merger
   triggers a burst.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="massRatio"   />
   <entry label="bulgeToTotal"/>
  </enumeration>
  !!]

  !![
  <mergerMassMovements name="mergerMassMovementsBaugh2005" docformat="rst">
   <description>
   A merger mass movements class which implements mass movements according to:

   * If :math:`M_\mathrm{satellite} &gt; f_\mathrm{major} M_\mathrm{central}` then all mass from both satellite and central galaxies moves to the spheroid :term:`component` of the central galaxy;
   * Otherwise:

     * If :math:`M_\mathrm{satellite} &gt; f_\mathrm{burst} M_\mathrm{central}` and the gas fraction in the host equals or exceeds :math:`f_\mathrm{gas,crit}` then all gas is moved to the host spheroid, while the host stellar disk remains in place. This is the criterion for a burst in a minor merger of :cite:t:`baugh_can_2005`. If ``[burstCriterion]``\ :math:`=`\ ``bulgeToTotal`` the first condition is instead :math:`M_\mathrm{central, spheroid} &lt; f_\mathrm{burst} M_\mathrm{central}`, as in earlier versions of Galacticus (``[burstCriterion]``\ :math:`=`\ ``massRatio``, the default, selects the criterion of :cite:t:`baugh_can_2005`).
     * Otherwise, gas from the satellite moves to the :term:`component` of the central specified by the ``[destinationGasMinorMerger]`` parameter (either "``disk``" or "``spheroid``"), stars from the satellite moves to the spheroid of the central and mass in the central does not move.

   Here, :math:`f_\mathrm{major}=`\ ``[massRatioMajorMerger]`` is the mass ratio above which a merger is considered to be "major", while :math:`f_\mathrm{burst}=`\ ``[ratioMassBurst]`` and :math:`f_\mathrm{gas,crit}=`\ ``[fractionGasCriticalBurst]``.
   </description>
  </mergerMassMovements>
  !!]
  type, extends(mergerMassMovementsMemoized) :: mergerMassMovementsBaugh2005
     !!{RST
     A merger mass movements class which uses the :cite:t:`baugh_can_2005` calculation.
     !!}
     private
     double precision                                         :: massRatioMajorMerger     , ratioMassBurst, &
          &                                                      fractionGasCriticalBurst
     type            (enumerationDestinationMergerType      ) :: destinationGasMinorMerger
     type            (enumerationBaugh2005BurstCriterionType) :: burstCriterion
   contains
     final     ::              baugh2005Destructor
     procedure :: calculate => baugh2005Calculate
  end type mergerMassMovementsBaugh2005

  interface mergerMassMovementsBaugh2005
     !!{RST
     Constructors for the :galacticus-class:`mergerMassMovementsBaugh2005` merger mass movements class.
     !!}
     module procedure baugh2005ConstructorParameters
     module procedure baugh2005ConstructorInternal
  end interface mergerMassMovementsBaugh2005

contains

  function baugh2005ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`mergerMassMovementsBaugh2005` merger mass movements class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerMassMovementsBaugh2005)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: massRatioMajorMerger     , ratioMassBurst, &
         &                                                           fractionGasCriticalBurst
    type            (varying_string              )                :: destinationGasMinorMerger, burstCriterion

    !![
    <inputParameter docformat="rst">
      <name>massRatioMajorMerger</name>
      <defaultValue>0.25d0</defaultValue>
      <description>
      The mass ratio above which mergers are considered to be "major". (Note that this default differs from the value,
      :math:`f_\mathrm{ellip}=0.3`, used by :cite:t:`baugh_can_2005`.)
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>ratioMassBurst</name>
      <defaultValue>0.05d0</defaultValue>
      <description>
      If ``[burstCriterion]`` is ``massRatio``, the satellite to central galaxy mass ratio above which minor mergers trigger a
      burst; if ``bulgeToTotal``, the bulge-to-total mass ratio of the central galaxy below which minor mergers trigger a burst.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>burstCriterion</name>
      <defaultValue>var_str('massRatio')</defaultValue>
      <description>
      The criterion used to decide whether a minor merger triggers a burst. If ``massRatio`` (the criterion of
      :cite:t:`baugh_can_2005`), a burst requires that the satellite to central galaxy mass ratio exceed ``[ratioMassBurst]``. If
      ``bulgeToTotal`` (the criterion used by earlier versions of Galacticus), a burst instead requires that the bulge-to-total
      mass ratio of the central galaxy be less than ``[ratioMassBurst]``. In both cases the gas fraction of the central galaxy must
      also equal or exceed ``[fractionGasCriticalBurst]``.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>fractionGasCriticalBurst</name>
      <defaultValue>0.75d0</defaultValue>
      <description>
      The host gas fraction above which mergers are considered to trigger a burst.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>destinationGasMinorMerger</name>
      <defaultValue>var_str('spheroid')</defaultValue>
      <description>
      The component to which satellite galaxy gas moves to as a result of a minor merger.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerMassMovementsBaugh2005(massRatioMajorMerger,enumerationDestinationMergerEncode(char(destinationGasMinorMerger),includesPrefix=.false.),ratioMassBurst,fractionGasCriticalBurst,enumerationBaugh2005BurstCriterionEncode(char(burstCriterion),includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function baugh2005ConstructorParameters

  function baugh2005ConstructorInternal(massRatioMajorMerger,destinationGasMinorMerger,ratioMassBurst,fractionGasCriticalBurst,burstCriterion) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`mergerMassMovementsBaugh2005` merger mass movements.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (mergerMassMovementsBaugh2005          )                :: self
    double precision                                        , intent(in   ) :: massRatioMajorMerger     , ratioMassBurst, &
         &                                                                     fractionGasCriticalBurst
    type            (enumerationDestinationMergerType      ), intent(in   ) :: destinationGasMinorMerger
    type            (enumerationBaugh2005BurstCriterionType), intent(in   ) :: burstCriterion
    !![
    <constructorAssign variables="massRatioMajorMerger, destinationGasMinorMerger, ratioMassBurst, fractionGasCriticalBurst, burstCriterion"/>
    !!]

    if (.not.enumerationBaugh2005BurstCriterionIsValid(burstCriterion)) call Error_Report('[burstCriterion] is invalid'//{introspection:location})

    return
  end function baugh2005ConstructorInternal

  subroutine baugh2005Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`mergerMassMovementsBaugh2005` merger mass movements class.
    !!}
    implicit none
    type(mergerMassMovementsBaugh2005), intent(inout) :: self

    call self%detachHooks()
    return
  end subroutine baugh2005Destructor

  subroutine baugh2005Calculate(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !!{RST
    Determine how different mass components should be redistributed as the result of a merger according to the model of :cite:t:`baugh_can_2005`.
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
    mergerIsMajor                   =    massSatellite    >= self%massRatioMajorMerger    *massHost
    select case (self%burstCriterion%ID)
    case (baugh2005BurstCriterionMassRatio   %ID)
       ! The criterion of Baugh et al. (2005): the satellite to central mass ratio must exceed the burst threshold.
       triggersBurst                =    massSatellite    >  self%ratioMassBurst          *massHost
    case (baugh2005BurstCriterionBulgeToTotal%ID)
       ! The criterion used by earlier versions of Galacticus: the central bulge-to-total ratio must be below the burst threshold.
       triggersBurst                =    massSpheroidHost <  self%ratioMassBurst          *massHost
    case default
       triggersBurst                =   .false.
    end select
    triggersBurst                   =    mergerIsMajor                                              &
         &                             .or.                                                         &
         &                              (                                                           &
         &                               triggersBurst                                              &
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
  end subroutine baugh2005Calculate
