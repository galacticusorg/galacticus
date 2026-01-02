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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorAgesStellarMassWeighted">
   <description>
     A node property extractor which extracts stellar mass-weighted ages for disk and spheroid components. Requires the
     \refClass{nodeOperatorAgesStellarMassWeighted} node operator to be used to accumulate the relevant integrals for each disk
     and spheroid.

     Specifically, the quantities computed by this class are     
     \begin{eqnarray}
     \langle t \rangle &amp;=&amp; \left. \int_0^t \mathrm{d}t^\prime (t-t^\prime) \dot{\psi}(t^\prime) \right. \int_0^t \mathrm{d}t^\prime \dot{\psi}(t^\prime) \nonumber \\
                       &amp;=&amp; t - \left. \int_0^t \mathrm{d}t^\prime t^\prime \dot{\psi}(t^\prime) \right. \int_0^t \mathrm{d}t^\prime \dot{\psi}(t^\prime),
     \end{eqnarray}
     where $\dot{\psi}(t^\prime)$ is the star formation rate at time $t^\prime$ and $t$ is the present time..
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorAgesStellarMassWeighted
     !!{
     A property extractor which extracts stellar mass-weighted ages of disk and spheroid.
     !!}
     private
     integer:: stellarMassFormedDiskID    , timeStellarMassFormedDiskID    , &
          &    stellarMassFormedSpheroidID, timeStellarMassFormedSpheroidID, &
          &    stellarMassFormedNSCID     , timeStellarMassFormedNSCID     
   contains
     procedure :: elementCount => agesStellarMassWeightedElementCount
     procedure :: extract      => agesStellarMassWeightedExtract
     procedure :: names        => agesStellarMassWeightedNames
     procedure :: descriptions => agesStellarMassWeightedDescriptions
     procedure :: unitsInSI    => agesStellarMassWeightedUnitsInSI
  end type nodePropertyExtractorAgesStellarMassWeighted

  interface nodePropertyExtractorAgesStellarMassWeighted
     !!{
     Constructors for the \refClass{nodePropertyExtractorAgesStellarMassWeighted} output extractor class.
     !!}
     module procedure agesStellarMassWeightedConstructorParameters
     module procedure agesStellarMassWeightedConstructorInternal
  end interface nodePropertyExtractorAgesStellarMassWeighted

contains

  function agesStellarMassWeightedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorAgesStellarMassWeighted} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorAgesStellarMassWeighted)                :: self
    type (inputParameters                             ), intent(inout) :: parameters

    self=nodePropertyExtractorAgesStellarMassWeighted()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function agesStellarMassWeightedConstructorParameters

  function agesStellarMassWeightedConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorAgesStellarMassWeighted} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorAgesStellarMassWeighted) :: self
    
    !![
    <addMetaProperty component="disk"     name="agesStellarMassFormed"     id="self%stellarMassFormedDiskID"         isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="disk"     name="agesTimeStellarMassFormed" id="self%timeStellarMassFormedDiskID"     isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="spheroid" name="agesStellarMassFormed"     id="self%stellarMassFormedSpheroidID"     isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="spheroid" name="agesTimeStellarMassFormed" id="self%timeStellarMassFormedSpheroidID" isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="NSC"      name="agesStellarMassFormed"     id="self%stellarMassFormedNSCID"          isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="NSC"      name="agesTimeStellarMassFormed" id="self%timeStellarMassFormedNSCID"      isEvolvable="yes" isCreator="no" />
    !!]
    return
  end function agesStellarMassWeightedConstructorInternal

  integer function agesStellarMassWeightedElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily agesStellarMassWeighted} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorAgesStellarMassWeighted), intent(inout) :: self
    double precision                                              , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    agesStellarMassWeightedElementCount=3
    return
  end function agesStellarMassWeightedElementCount

  function agesStellarMassWeightedExtract(self,node,time,instance)
    !!{
    Implement a agesStellarMassWeighted output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    implicit none
    double precision                                              , dimension(:) , allocatable :: agesStellarMassWeightedExtract
    class           (nodePropertyExtractorAgesStellarMassWeighted), intent(inout), target      :: self
    type            (treeNode                                    ), intent(inout), target      :: node
    double precision                                              , intent(in   )              :: time
    type            (multiCounter                                ), intent(inout), optional    :: instance
    class           (nodeComponentDisk                           )               , pointer     :: disk
    class           (nodeComponentSpheroid                       )               , pointer     :: spheroid
    class           (nodeComponentNSC                            )               , pointer     :: nuclearStarCluster

    double precision                                                                           :: massStellarDisk               , massTimeStellarDisk              , &
         &                                                                                        massStellarSpheroid           , massTimeStellarSpheroid          , &
         &                                                                                        massStellarNuclearStarCluster , massTimeStellarNuclearStarCluster, &
         &                                                                                        ageDisk                       , ageSpheroid                      , &
         &                                                                                        ageNuclearStarCluster              
    !$GLC attributes unused :: time, instance

    ! Extract required quantities.
    disk               => node%disk    ()
    spheroid           => node%spheroid()
    nuclearStarCluster => node%NSC     ()

    select type (disk              )
    type is (nodeComponentDisk    )
       ! Disk does not yet exist.
       massStellarDisk                   =0.0d0
       massTimeStellarDisk               =0.0d0
    class default
       massStellarDisk                   =disk               %floatRank0MetaPropertyGet(self%    stellarMassFormedDiskID    )
       massTimeStellarDisk               =disk               %floatRank0MetaPropertyGet(self%timeStellarMassFormedDiskID    )
    end select

    select type (spheroid          )
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist.
       massStellarSpheroid               =0.0d0
       massTimeStellarSpheroid           =0.0d0
    class default
       massStellarSpheroid               =spheroid           %floatRank0MetaPropertyGet(self%    stellarMassFormedSpheroidID)
       massTimeStellarSpheroid           =spheroid           %floatRank0MetaPropertyGet(self%timeStellarMassFormedSpheroidID)
    end select

    select type (nuclearStarCluster)
    type is (nodeComponentNSC     )
        ! Nuclear star cluster does not yet exist.
        massStellarNuclearStarCluster    =0.0d0
        massTimeStellarNuclearStarCluster=0.0d0
    class default
        massStellarNuclearStarCluster    =nuclearStarCluster%floatRank0MetaPropertyGet(self%    stellarMassFormedNSCID      )
        massTimeStellarNuclearStarCluster=nuclearStarCluster%floatRank0MetaPropertyGet(self%timeStellarMassFormedNSCID      ) 
    end select     
    ! Compute ages.

    if (massStellarDisk > 0.0d0) then
       ageDisk              =+time                              &
            &                -massTimeStellarDisk               &
            &                /massStellarDisk
    else
       ageDisk              =-1.0d0
    end if

    if (massStellarSpheroid > 0.0d0) then
       ageSpheroid          =+time                              &
            &                -massTimeStellarSpheroid           &
            &                /massStellarSpheroid
    else
       ageSpheroid          =-1.0d0
    end if

    if (massStellarNuclearStarCluster > 0.0d0) then
       ageNuclearStarCluster=+time                              &
            &                -massTimeStellarNuclearStarCluster &
            &                /massStellarNuclearStarCluster
    else
       ageNuclearStarCluster     =-1.0d0
    end if 

    ! Set return results.
    allocate(agesStellarMassWeightedExtract(3))
    agesStellarMassWeightedExtract=[                       &
         &                          ageDisk              , &
         &                          ageSpheroid          , &
         &                          ageNuclearStarCluster  &
         &                         ]
    return
  end function agesStellarMassWeightedExtract

  subroutine agesStellarMassWeightedNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily agesStellarMassWeighted} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorAgesStellarMassWeighted), intent(inout)                             :: self
    double precision                                              , intent(in   )                             :: time
    type            (varying_string                              ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(3))
    names(1)=var_str('diskAgeStellarMassWeighted'              )
    names(2)=var_str('spheroidAgeStellarMassWeighted'          )
    names(3)=var_str('nuclearStarClusterAgeStellarMassWeighted')
    return
  end subroutine agesStellarMassWeightedNames

  subroutine agesStellarMassWeightedDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily agesStellarMassWeighted} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorAgesStellarMassWeighted), intent(inout)                             :: self
    double precision                                              , intent(in   )                             :: time
    type            (varying_string                              ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)=var_str('Stellar mass-weighted age of the disk [Gyr].'                )
    descriptions(2)=var_str('Stellar mass-weighted age of the spheroid [Gyr].'            )
    descriptions(3)=var_str('stellar mass-weighted age of the nuclear star cluster [Gyr].')
    return
  end subroutine agesStellarMassWeightedDescriptions

  function agesStellarMassWeightedUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily agesStellarMassWeighted} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    double precision                                              , dimension(:) , allocatable :: agesStellarMassWeightedUnitsInSI
    class           (nodePropertyExtractorAgesStellarMassWeighted), intent(inout)              :: self
    double precision                                              , intent(in   )              :: time
   !$GLC attributes unused :: self, time

    allocate(agesStellarMassWeightedUnitsInSI(3))
    agesStellarMassWeightedUnitsInSI=gigaYear
    return
  end function agesStellarMassWeightedUnitsInSI

