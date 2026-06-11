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

  !!{RST
  Implements a property extractor class for the core mass of the :term:`FDM` soliton.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorSoliton" docformat="rst">
   <description>
   Extracts physical properties of the fuzzy dark matter (:term:`FDM`) soliton core (such as core radius and core mass) associated with each halo node, enabling analysis of quantum pressure effects in FDM models.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorSoliton
     !!{RST
     A property extractor class for the properties of the :term:`FDM` soliton.
     !!}
     private
     integer :: massCoreNormalID, massCoreID, densityCoreID, radiusCoreID, radiusSolitonID, zetaID
   contains
     procedure :: elementCount => solitonElementCount
     procedure :: extract      => solitonExtract
     procedure :: names        => solitonNames
     procedure :: descriptions => solitonDescriptions
     procedure :: unitsInSI    => solitonUnitsInSI
     procedure :: units       => solitonUnits
  end type nodePropertyExtractorSoliton

  interface nodePropertyExtractorSoliton
     !!{RST
     Constructors for the ``nodePropertyExtractorSoliton`` property extractor class.
     !!}
     module procedure solitonConstructorParameters
     module procedure solitonConstructorInternal
  end interface nodePropertyExtractorSoliton

contains

  function solitonConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``nodePropertyExtractorSoliton`` property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorSoliton)                :: self
    type(inputParameters             ), intent(inout) :: parameters

    self=nodePropertyExtractorSoliton()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function solitonConstructorParameters

  function solitonConstructorInternal() result(self)
    !!{RST
    Internal constructor for the ``nodePropertyExtractorSoliton`` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSoliton) :: self
    
    !![
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreNormal" id="self%massCoreNormalID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"       id="self%massCoreID"       isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonDensityCore"    id="self%densityCoreID"    isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusCore"     id="self%radiusCoreID"     isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusSoliton"  id="self%radiusSolitonID"  isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonZeta"           id="self%zetaID"           isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function solitonConstructorInternal

  integer function solitonElementCount(self,time)
    !!{RST
    Return the number of elements in the ``soliton`` property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorSoliton), intent(inout) :: self
    double precision                              , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    solitonElementCount=6
    return
  end function solitonElementCount

  function solitonExtract(self,node,time,instance)
    !!{RST
    Implement a ``soliton`` property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    double precision                     , dimension(: ), allocatable :: solitonExtract
    class(nodePropertyExtractorSoliton  ), intent(inout), target      :: self
    type (treeNode                      ), intent(inout), target      :: node
    double precision                     , intent(in   )              :: time
    type (multiCounter                  ), intent(inout), optional    :: instance
    class(nodeComponentDarkMatterProfile)               , pointer     :: darkMatterProfile
    !$GLC attributes unused :: time, instance

    allocate(solitonExtract(6))
    darkMatterProfile => node%darkMatterProfile()
    select type (darkMatterProfile)
    type is (nodeComponentDarkMatterProfile)
       ! Dark matter profile does not exist.
      solitonExtract=[        & 
        &             0.0d0,  &
        &             0.0d0,  &
        &             0.0d0,  &
        &             0.0d0,  &
        &             0.0d0,  &
        &             0.0d0   &
        &            ]
    class default
      solitonExtract=[                                                                    &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreNormalID), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreID      ), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%densityCoreID   ), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%radiusCoreID    ), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%radiusSolitonID ), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%zetaID          )  &
       &             ]    
    end select
    return
  end function solitonExtract

  subroutine solitonNames(self,time,names)
    !!{RST
    Return the names of the ``soliton`` property.
    !!}
    implicit none
    class(nodePropertyExtractorSoliton), intent(inout)                             :: self
    double precision                   , intent(in   )                             :: time
    type            (varying_string   ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(6))
    names(1)=var_str('solitonMassCoreNormal')
    names(2)=var_str('solitonMassCore'      )
    names(3)=var_str('solitonDensityCore'   )
    names(4)=var_str('solitonRadiusCore'    )
    names(5)=var_str('solitonRadiusSoliton' )
    names(6)=var_str('solitonZetazOverZeta0')
    return
  end subroutine solitonNames

  subroutine solitonDescriptions(self,time,descriptions)
    !!{RST
    Return the descriptions of the ``soliton`` property.
    !!}
    implicit none
    class(nodePropertyExtractorSoliton), intent(inout)                             :: self
    double precision                   , intent(in   )                             :: time
    type            (varying_string   ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(6))
    descriptions(1)=var_str('The solitonic core mass of the FDM halo (without scatter), in units of M☉.')
    descriptions(2)=var_str('The solitonic core mass of the FDM halo, in units of M☉.'                  )
    descriptions(3)=var_str('The soliton central density of the FDM halo, in units of [M☉/Mpc³].'       )
    descriptions(4)=var_str('The soliton core radius of the FDM halo, in units of Mpc.'                 )
    descriptions(5)=var_str('The soliton transition radius of the FDM halo, in units of Mpc.'           )
    descriptions(6)=var_str('The ratio of density contrast at redshift z and 0, dimensionless.'         )
    return
  end subroutine solitonDescriptions

  function solitonUnitsInSI(self,time)
    !!{RST
    Return the units of the ``Soliton`` property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                   , allocatable  , dimension(:) :: solitonUnitsInSI
    class(nodePropertyExtractorSoliton), intent(inout)               :: self
    double precision                                 , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(solitonUnitsInSI(6))
    solitonUnitsInSI=[                         &
         &            massSolar              , &
         &            massSolar              , &
         &            massSolar/megaParsec**3, &
         &            megaParsec             , &
         &            megaParsec             , &
         &            1.0d0                    &
         &           ]
    return
  end function solitonUnitsInSI

  function solitonUnits(self,time) result(units)
    !!{RST
    Return the units of the soliton properties.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType                    ), dimension(:) , allocatable :: units
    class           (nodePropertyExtractorSoliton), intent(inout)              :: self
    double precision                              , intent(in   )              :: time
    double precision                              , dimension(:) , allocatable :: siValues
    !$GLC attributes unused :: self

    siValues=self%unitsInSI(time)
    allocate(units(6))
    units(1)=unitType(siValues(1),description='Solar masses',quantity='solMass'      )
    units(2)=unitType(siValues(2),description='Solar masses',quantity='solMass'      )
    units(3)=unitType(siValues(3),description='M☉/Mpc³'     ,quantity='solMass/Mpc^3')
    units(4)=unitType(siValues(4),description='Mpc'         ,quantity='Mpc'          )
    units(5)=unitType(siValues(5),description='Mpc'         ,quantity='Mpc'          )
    units(6)=unitType(siValues(6)                                                    )
    return
  end function solitonUnits
