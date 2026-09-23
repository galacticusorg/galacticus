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
  Implements a merger mass movements class which uses a simple calculation.
  !!}

  !![
  <mergerMassMovements name="mergerMassMovementsSimple" docformat="rst">
   <description>
   A merger mass movements class which implements mass movements according to:

   * If :math:`\min(M_\mathrm{satellite},M_\mathrm{central}) \ge f_\mathrm{major} \max(M_\mathrm{satellite},M_\mathrm{central})` then all mass from both satellite and central galaxies moves to the spheroid :term:`component` of the central galaxy;
   * Otherwise: gas and stars from the satellite move to the :term:`component`\ s of the central specified by the ``[destinationGasMinorMerger]`` and ``[destinationStarsMinorMerger]`` parameters respectively (each either "``disk``", "``spheroid``", or "``dominant``"---the latter meaning whichever of the disk and spheroid of the more massive galaxy is the more massive), and mass in the central does not move.

   Here, :math:`f_\mathrm{major}=`\ ``[massRatioMajorMerger]`` is the mass ratio above which a merger is considered to be "major". Note that the masses used in this criterion are the total galactic masses (gas plus stars, in both disk and spheroid) of the two galaxies.
   </description>
  </mergerMassMovements>
  !!]
  type, extends(mergerMassMovementsMemoized) :: mergerMassMovementsSimple
     !!{RST
     A merger mass movements class which uses a simple calculation.
     !!}
     private
     double precision                                   :: massRatioMajorMerger
     type            (enumerationDestinationMergerType) :: destinationGasMinorMerger, destinationStarsMinorMerger
   contains
     final     ::              simpleDestructor
     procedure :: calculate => simpleCalculate
  end type mergerMassMovementsSimple

  interface mergerMassMovementsSimple
     !!{RST
     Constructors for the :galacticus-class:`mergerMassMovementsSimple` merger mass movements class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface mergerMassMovementsSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`mergerMassMovementsSimple` merger mass movements class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerMassMovementsSimple)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: massRatioMajorMerger
    type            (varying_string           )                :: destinationGasMinorMerger, destinationStarsMinorMerger

    !![
    <inputParameter docformat="rst">
      <name>massRatioMajorMerger</name>
      <defaultValue>0.25d0</defaultValue>
      <description>
      The mass ratio above which mergers are considered to be "major".
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
    <inputParameter docformat="rst">
      <name>destinationStarsMinorMerger</name>
      <defaultValue>var_str('spheroid')</defaultValue>
      <description>
      The component to which satellite galaxy stars move to as a result of a minor merger.
      </description>
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
    !!{RST
    Internal constructor for the :galacticus-class:`mergerMassMovementsSimple` merger mass movements class.
    !!}
    implicit none
    type            (mergerMassMovementsSimple       )                        :: self
    double precision                                  , intent(in   )         :: massRatioMajorMerger
    type            (enumerationDestinationMergerType), intent(in   )         :: destinationGasMinorMerger, destinationStarsMinorMerger
    !![
    <constructorAssign variables="massRatioMajorMerger, destinationGasMinorMerger, destinationStarsMinorMerger"/>
    !!]

    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`mergerMassMovementsSimple` merger mass movements class
    !!}
    implicit none
    type(mergerMassMovementsSimple), intent(inout) :: self

    call self%detachHooks()
    return
  end subroutine simpleDestructor

  subroutine simpleCalculate(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !!{RST
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

    nodeHost                 => node                     %mergesWith      (                         )
    massDistributionHost     => nodeHost                 %massDistribution(massType=massTypeGalactic)
    massDistributionSatellite=> node                     %massDistribution(massType=massTypeGalactic)
    massSatellite            =  massDistributionSatellite%massTotal       (                         )
    massHost                 =  massDistributionHost     %massTotal       (                         )
    mergerIsMajor            =  massSatellite > 0.0d0 .and. massHost > 0.0d0 .and. min(massSatellite,massHost) >= self%massRatioMajorMerger*max(massSatellite,massHost)
    !![
    <objectDestructor name="massDistributionHost"     />
    <objectDestructor name="massDistributionSatellite"/>
    !!]
    if (mergerIsMajor) then
       destinationGasSatellite  =     destinationMergerSpheroid
       destinationStarsSatellite=     destinationMergerSpheroid
       destinationGasHost       =     destinationMergerSpheroid
       destinationStarsHost     =     destinationMergerSpheroid
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
          destinationGasSatellite=     destinationDominant
       else
          destinationGasSatellite=self%destinationGasMinorMerger
       end if
       if (self%destinationStarsMinorMerger == destinationMergerDominant) then
          destinationStarsSatellite=     destinationDominant
       else
          destinationStarsSatellite=self%destinationStarsMinorMerger
       end if
       destinationGasHost  =     destinationMergerUnmoved
       destinationStarsHost=     destinationMergerUnmoved
    end if
    return
  end subroutine simpleCalculate
