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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorStarFormationRateInterOutput">
   <description>
    A node property extractor which extracts the mean star formation rate between successive outputs. Intended to be paired with the
    \refClass{nodeOperatorStarFormationRateInterOutput} class to compute those rates.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorStarFormationRateInterOutput
     !!{
     A property extractor which extracts mean star formation rates between successive outputs.
     !!}
     private
     integer :: starFormationRateDiskInterOutputID, starFormationRateSpheroidInterOutputID, &
          &     starFormationRateNSCInterOutputID , starFormationRateInterOutputNextID
   contains
     procedure :: elementCount => starFormationRateInterOutputElementCount
     procedure :: extract      => starFormationRateInterOutputExtract
     procedure :: names        => starFormationRateInterOutputNames
     procedure :: descriptions => starFormationRateInterOutputDescriptions
     procedure :: unitsInSI    => starFormationRateInterOutputUnitsInSI
  end type nodePropertyExtractorStarFormationRateInterOutput

  interface nodePropertyExtractorStarFormationRateInterOutput
     !!{
     Constructors for the \refClass{nodePropertyExtractorStarFormationRateInterOutput} output extractor class.
     !!}
     module procedure starFormationRateInterOutputConstructorParameters
     module procedure starFormationRateInterOutputConstructorInternal
  end interface nodePropertyExtractorStarFormationRateInterOutput

contains

  function starFormationRateInterOutputConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorStarFormationRateInterOutput} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorStarFormationRateInterOutput)                :: self
    type(inputParameters                                  ), intent(inout) :: parameters

    self=nodePropertyExtractorStarFormationRateInterOutput()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function starFormationRateInterOutputConstructorParameters

  function starFormationRateInterOutputConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorStarFormationRateInterOutput} output extractor property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorStarFormationRateInterOutput) :: self
    
    !![
    <addMetaProperty component="disk"     name="starFormationRateDiskInterOutput" isEvolvable="yes" id="self%starFormationRateDiskInterOutputID"    />
    <addMetaProperty component="spheroid" name="starFormationRateDiskInterOutput" isEvolvable="yes" id="self%starFormationRateSpheroidInterOutputID"/>
    <addMetaProperty component="NSC"      name="starFormationRateNSCInterOutput"  isEvolvable="yes" id="self%starFormationRateNSCInterOutputID"     />
    !!]
    return
  end function starFormationRateInterOutputConstructorInternal

  integer function starFormationRateInterOutputElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily starFormationRateInterOutput} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorStarFormationRateInterOutput), intent(inout) :: self
    double precision                                                   , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    starFormationRateInterOutputElementCount=3
    return
  end function starFormationRateInterOutputElementCount

  function starFormationRateInterOutputExtract(self,node,time,instance)
    !!{
    Implement a starFormationRateInterOutput output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    implicit none
    double precision                                                   , dimension(:) , allocatable :: starFormationRateInterOutputExtract
    class           (nodePropertyExtractorStarFormationRateInterOutput), intent(inout), target      :: self
    type            (treeNode                                         ), intent(inout), target      :: node
    double precision                                                   , intent(in   )              :: time
    type            (multiCounter                                     ), intent(inout), optional    :: instance
    class           (nodeComponentDisk                                ), pointer                    :: disk
    class           (nodeComponentSpheroid                            ), pointer                    :: spheroid
    class           (nodeComponentNSC                                 ), pointer                    :: nuclearStarCluster
    double precision                                                                                :: starFormationRateDisk              , starFormationRateSpheroid, &
                   &                                                                                   starFormationRatenuclearStarCluster
    !$GLC attributes unused :: time, instance

    allocate(starFormationRateInterOutputExtract(3))
    disk               => node%disk    ()
    spheroid           => node%spheroid()
    nuclearStarCluster => node%NSC     ()
    select type (disk    )
    type is (nodeComponentDisk    )
       ! Disk does not yet exist.
       starFormationRateDisk               =0.0d0
    class default
       starFormationRateDisk               =disk               %floatRank0MetaPropertyGet(self%starFormationRateDiskInterOutputID    )
    end select
    select type (spheroid)
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist.
       starFormationRateSpheroid           =0.0d0
    class default
       starFormationRateSpheroid           =spheroid           %floatRank0MetaPropertyGet(self%starFormationRateSpheroidInterOutputID)
    end select
    select type (nuclearStarCluster)
    type is (nodeComponentNSC     )
      ! NSC does not yet exist.
        starFormationRatenuclearStarCluster=0.0d0
    class default
        starFormationRatenuclearStarCluster=nuclearStarCluster%floatRank0MetaPropertyGet(self%starFormationRateNSCInterOutputID     )
    end select

    starFormationRateInterOutputExtract=[                                     &
         &                               starFormationRateDisk              , &
         &                               starFormationRateSpheroid          , &
         &                               starFormationRatenuclearStarCluster  &
         &                              ]
    return
  end function starFormationRateInterOutputExtract

  subroutine starFormationRateInterOutputNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily starFormationRateInterOutput} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorStarFormationRateInterOutput), intent(inout)                            :: self
    double precision                                                   , intent(in   )                            :: time
    type            (varying_string                                   ), intent(inout), dimension(:), allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(3))
    names(1)=var_str(              'diskStarFormationRateInterOutputMean')
    names(2)=var_str(          'spheroidStarFormationRateInterOutputMean')
    names(3)=var_str('nuclearStarClusterStarFormationRateInterOutputMean')
    return
  end subroutine starFormationRateInterOutputNames

  subroutine starFormationRateInterOutputDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily starFormationRateInterOutput} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorStarFormationRateInterOutput), intent(inout)                             :: self
    double precision                                                   , intent(in   )                             :: time
    type            (varying_string                                   ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)=var_str('Mean star formation rate in the disk between this output and the previous output.'                )
    descriptions(2)=var_str('Mean star formation rate in the spheroid between this output and the previous output.'            )
    descriptions(3)=var_str('Mean star formation rate in the nuclear star cluster between this output and the previous output.')
    return
  end subroutine starFormationRateInterOutputDescriptions

  function starFormationRateInterOutputUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily starFormationRateInterOutput} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, massSolar
    implicit none
    double precision                                                   , dimension(:) , allocatable :: starFormationRateInterOutputUnitsInSI
    class           (nodePropertyExtractorStarFormationRateInterOutput), intent(inout)              :: self
    double precision                                                   , intent(in   )              :: time
   !$GLC attributes unused :: self, time

    allocate(starFormationRateInterOutputUnitsInSI(3))
    starFormationRateInterOutputUnitsInSI=+massSolar &
         &                                /gigaYear
    return
  end function starFormationRateInterOutputUnitsInSI

