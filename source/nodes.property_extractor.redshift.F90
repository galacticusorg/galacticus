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
Implements a redshift property extractor class.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRedshift" docformat="rst">
   <description>
   Extracts the current cosmological redshift at which a node exists, computed from the node's cosmic time, and outputs it as the scalar property named "``redshift``".
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRedshift
     !!{RST
     A redshift property extractor class.
     !!}
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
   contains
     final     ::                redshiftDestructor
     procedure :: extract     => redshiftExtract
     procedure :: name        => redshiftName
     procedure :: description => redshiftDescription
     procedure :: unitsInSI   => redshiftUnitsInSI
     procedure :: units       => redshiftUnits
  end type nodePropertyExtractorRedshift

  interface nodePropertyExtractorRedshift
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorRedshift` property extractor class.
     !!}
     module procedure redshiftConstructorParameters
     module procedure redshiftConstructorInternal
  end interface nodePropertyExtractorRedshift

contains

  function redshiftConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorRedshift` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRedshift)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRedshift(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function redshiftConstructorParameters

  function redshiftConstructorInternal(cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorRedshift` property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRedshift)                        :: self
    class(cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    return
  end function redshiftConstructorInternal

  subroutine redshiftDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorRedshift` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRedshift), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine redshiftDestructor

  double precision function redshiftExtract(self,node,instance)
    !!{RST
    Implement a last isolated redshift output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nodePropertyExtractorRedshift), intent(inout), target   :: self
    type (treeNode                     ), intent(inout), target   :: node
    type (multiCounter                 ), intent(inout), optional :: instance
    class(nodeComponentBasic           ), pointer                 :: basic
    !$GLC attributes unused :: self, instance

    basic           => node %basic                                            ()
    redshiftExtract =  self %cosmologyFunctions_%redshiftFromExpansionFactor(    &
         &                         self %cosmologyFunctions_%expansionFactor (   &
         &                         basic                    %time             () &
         &                                                                   )   &
         &                                                                  )
    return
  end function redshiftExtract

  function redshiftName(self)
    !!{RST
    Return the name of the last isolated redshift property.
    !!}
    implicit none
    type (varying_string               )                :: redshiftName
    class(nodePropertyExtractorRedshift), intent(inout) :: self
    !$GLC attributes unused :: self

    redshiftName=var_str('redshift')
    return
  end function redshiftName

  function redshiftDescription(self)
    !!{RST
    Return a description of the redshift property.
    !!}
    implicit none
    type (varying_string               )                :: redshiftDescription
    class(nodePropertyExtractorRedshift), intent(inout) :: self
    !$GLC attributes unused :: self

    redshiftDescription=var_str('The redshift of the epoch at which the galaxy exists.')
    return
  end function redshiftDescription

  double precision function redshiftUnitsInSI(self)
    !!{RST
    Return the units of the last isolated redshift property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorRedshift), intent(inout) :: self
    !$GLC attributes unused :: self

    redshiftUnitsInSI=1.0d0
    return
  end function redshiftUnitsInSI

  function redshiftUnits(self) result(units)
    !!{RST
    Return the units of the redshift property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                     )                :: units
    class(nodePropertyExtractorRedshift), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(1.0d0)
    return
  end function redshiftUnits
