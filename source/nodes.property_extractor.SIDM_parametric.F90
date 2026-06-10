!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  !+    Contributions to this file made by: Niusha Ahvazi

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSIDMParametric">
   <description>
     A node property extractor which extracts dark matter profile properties for the SIDM parametric model.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorSIDMParametric
     !!{
     A property extractor which extracts dark matter profile properties for the SIDM parametric model.
     !!}
     private
     integer:: tauID, VmaxSIDMID, RmaxSIDMID, RhosSIDMID, RsSIDMID, RcSIDMID
   contains
     procedure :: elementCount => SIDMParametricElementCount
     procedure :: extract      => SIDMParametricExtract
     procedure :: names        => SIDMParametricNames
     procedure :: descriptions => SIDMParametricDescriptions
     procedure :: unitsInSI    => SIDMParametricUnitsInSI
  end type nodePropertyExtractorSIDMParametric

  interface nodePropertyExtractorSIDMParametric
     !!{
     Constructors for the \refClass{nodePropertyExtractorSIDMParametric} class.
     !!}
     module procedure SIDMParametricConstructorParameters
     module procedure SIDMParametricConstructorInternal
  end interface nodePropertyExtractorSIDMParametric

contains

  function SIDMParametricConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorSIDMParametric} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorSIDMParametric)                :: self
    type (inputParameters                    ), intent(inout) :: parameters

    self=nodePropertyExtractorSIDMParametric()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function SIDMParametricConstructorParameters

  function SIDMParametricConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorSIDMParametric} class.
    !!}
    implicit none
    type(nodePropertyExtractorSIDMParametric) :: self
    
    !![
    <addMetaProperty component="darkMatterProfile" name="tau"      id="self%tauID"      isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="VmaxSIDM" id="self%VmaxSIDMID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="RmaxSIDM" id="self%RmaxSIDMID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="RhosSIDM" id="self%RhosSIDMID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="RsSIDM"   id="self%RsSIDMID"   isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="RcSIDM"   id="self%RcSIDMID"   isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function SIDMParametricConstructorInternal

  integer function SIDMParametricElementCount(self,time)
    !!{
    Return the number of elements in the \mono{SIDMParametric} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorSIDMParametric), intent(inout) :: self
    double precision                                     , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    SIDMParametricElementCount=6
    return
  end function SIDMParametricElementCount

  function SIDMParametricExtract(self,node,time,instance)
    !!{
    Extract parameters of the SIDM parametric model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkmatterProfile 
    implicit none
    double precision                                     , dimension(:) , allocatable :: SIDMParametricExtract
    class           (nodePropertyExtractorSIDMParametric), intent(inout), target      :: self
    type            (treeNode                           ), intent(inout), target      :: node
    double precision                                     , intent(in   )              :: time
    type            (multiCounter                       ), intent(inout), optional    :: instance
    class           (nodeComponentBasic                 )               , pointer     :: basic
    class           (nodeComponentDarkmatterProfile     )               , pointer     :: darkMatterProfile
    double precision                                                                  :: tauSIDMParametric     , VmaxSIDMParametric, &
         &                                                                               RmaxSIDMParametric    , RhosSIDMParametric, &
         &                                                                               RsSIDMParametric      , RcSIDMParametric
    !$GLC attributes unused :: time, instance

    ! Extract required quantities.
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    select type (darkMatterProfile)
    type is (nodeComponentDarkMatterProfile)
       ! The dark matter profile component does not yet exist, so return zero values.
       tauSIDMParametric =0.0d0
       VmaxSIDMParametric=0.0d0
       RmaxSIDMParametric=0.0d0  
       RhosSIDMParametric=0.0d0
       RsSIDMParametric  =0.0d0
       RcSIDMParametric  =0.0d0
    class default
       tauSIDMParametric =darkMatterProfile%floatRank0MetaPropertyGet(self%tauID     )
       VmaxSIDMParametric=darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID)
       RmaxSIDMParametric=darkMatterProfile%floatRank0MetaPropertyGet(self%RmaxSIDMID)
       RhosSIDMParametric=darkMatterProfile%floatRank0MetaPropertyGet(self%RhosSIDMID)
       RsSIDMParametric  =darkMatterProfile%floatRank0MetaPropertyGet(self%RsSIDMID  )
       RcSIDMParametric  =darkMatterProfile%floatRank0MetaPropertyGet(self%RcSIDMID  )
    end select
    ! Set return results.
    allocate(SIDMParametricExtract(6))
    SIDMParametricExtract=[                    &
         &                 tauSIDMParametric , &
         &                 VmaxSIDMParametric, &
         &                 RmaxSIDMParametric, &
         &                 RhosSIDMParametric, &
         &                 RsSIDMParametric  , &
         &                 RcSIDMParametric    &
         &                ]
    return
  end function SIDMParametricExtract

  subroutine SIDMParametricNames(self,time,names)
    !!{
    Return the names of the \mono{SIDMParametric} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorSIDMParametric), intent(inout)                             :: self
    double precision                                     , intent(in   )                             :: time
    type            (varying_string                     ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(6))
    names(1)=var_str('darkMatterProfileSIDMParametrictau' )
    names(2)=var_str('darkMatterProfileSIDMParametricVmax')
    names(3)=var_str('darkMatterProfileSIDMParametricRmax')
    names(4)=var_str('darkMatterProfileSIDMParametricRhos')
    names(5)=var_str('darkMatterProfileSIDMParametricRs'  )
    names(6)=var_str('darkMatterProfileSIDMParametricRc'  )
    return
  end subroutine SIDMParametricNames

  subroutine SIDMParametricDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the \mono{SIDMParametric} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorSIDMParametric), intent(inout)                            :: self
    double precision                                     , intent(in   )                            :: time
    type            (varying_string                     ), intent(inout), dimension(:), allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(6))
    descriptions(1)=var_str('Dimensionless time variable.'                                                                                     )
    descriptions(2)=var_str('Maximum velocity of the dark matter profile assuming the SIDM Parametric model [km/s].'                           )
    descriptions(3)=var_str('Radius at which the maximum velocity of the dark matter profile occurs, assuming the SIDM Parametric model [Mpc].')
    descriptions(4)=var_str('Scale density (SIDM) [M☉/Mpc³].'                                                                                  )
    descriptions(5)=var_str('Scale radius (SIDM) [Mpc].'                                                                                       )
    descriptions(6)=var_str('Core size [Mpc].'                                                                                                 )
    return
  end subroutine SIDMParametricDescriptions

  function SIDMParametricUnitsInSI(self,time)
    !!{
    Return the units of the \mono{SIDMParametric} properties in the SI system.
    !!}
    use :: Numerical_Constants_Prefixes,     only : kilo
    use :: Numerical_Constants_Astronomical, only : megaParsec, massSolar
    implicit none
    double precision                                     , dimension(:) , allocatable :: SIDMParametricUnitsInSI
    class           (nodePropertyExtractorSIDMParametric), intent(inout)              :: self
    double precision                                     , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(SIDMParametricUnitsInSI(6))
    SIDMParametricUnitsInSI(1)=1.0d0                   ! τ is dimensionless.
    SIDMParametricUnitsInSI(2)=kilo
    SIDMParametricUnitsInSI(3)=          megaParsec
    SIDMParametricUnitsInSI(4)=massSolar/megaParsec**3
    SIDMParametricUnitsInSI(5)=          megaParsec
    SIDMParametricUnitsInSI(6)=          megaParsec
    return
  end function SIDMParametricUnitsInSI

