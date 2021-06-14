!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements a halo mass output analysis property extractor class.

  use :: Virial_Density_Contrast, only : virialDensityContrastClass

  !# <nodePropertyExtractor name="nodePropertyExtractorMassHalo">
  !#  <description>A halo mass output analysis property extractor class.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassHalo
     !% A halo mass property extractor output analysis class. The property extracted is the ''\gls{dmou}'' mass of the halo within
     !% a radius enclosing a density contrast as defined by the supplied {\normalfont \ttfamily virialDensityContrast} class
     !% object. Note that the density contrast is defined here at the time at which the halo presently exists, \emph{not} at the
     !% time at which is was last isolated (as is used for standard definition of halo mass).
     private
     class(virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
   contains
     final     ::                massHaloDestructor
     procedure :: extract     => massHaloExtract
     procedure :: type        => massHaloType
     procedure :: name        => massHaloName
     procedure :: description => massHaloDescription
     procedure :: unitsInSI   => massHaloUnitsInSI
  end type nodePropertyExtractorMassHalo

  interface nodePropertyExtractorMassHalo
     !% Constructors for the ``massHalo'' output analysis class.
     module procedure massHaloConstructorParameters
     module procedure massHaloConstructorInternal
  end interface nodePropertyExtractorMassHalo

contains

  function massHaloConstructorParameters(parameters) result(self)
    !% Constructor for the ``massHalo'' output analysis property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorMassHalo)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(virialDensityContrastClass   ), pointer       :: virialDensityContrast_

    !# <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    self=nodePropertyExtractorMassHalo(virialDensityContrast_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="virialDensityContrast_"/>
    return
  end function massHaloConstructorParameters

  function massHaloConstructorInternal(virialDensityContrast_) result(self)
    !% Internal constructor for the ``massHalo'' output analysis property extractor class.
    implicit none
    type (nodePropertyExtractorMassHalo)                        :: self
    class(virialDensityContrastClass   ), intent(in   ), target :: virialDensityContrast_
    !# <constructorAssign variables="*virialDensityContrast_"/>

    return
  end function massHaloConstructorInternal

  subroutine massHaloDestructor(self)
    !% Destructor for the {\normalfont \ttfamily mass} output analysis property extractor class.
    implicit none
    type(nodePropertyExtractorMassHalo), intent(inout) :: self

    !# <objectDestructor name="self%virialDensityContrast_"/>
    return
  end subroutine massHaloDestructor

  double precision function massHaloExtract(self,node,instance)
    !% Implement a massHalo output analysis.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class(nodePropertyExtractorMassHalo), intent(inout)           :: self
    type (treeNode                     ), intent(inout), target   :: node
    type (multiCounter                 ), intent(inout), optional :: instance
    class(nodeComponentBasic           ), pointer                 :: basic
    !$GLC attributes unused :: instance

    basic           => node%basic()
    massHaloExtract =  Dark_Matter_Profile_Mass_Definition(node,self%virialDensityContrast_%densityContrast(basic%mass(),basic%time()))
    return
  end function massHaloExtract

  integer function massHaloType(self)
    !% Return the type of the halo mass property.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorMassHalo), intent(inout) :: self
    !$GLC attributes unused :: self

    massHaloType=outputAnalysisPropertyTypeLinear
    return
  end function massHaloType

  function massHaloName(self)
    !% Return the name of the massHalo property.
    implicit none
    type (varying_string               )                :: massHaloName
    class(nodePropertyExtractorMassHalo), intent(inout) :: self
    !$GLC attributes unused :: self

    massHaloName=var_str('massHaloEnclosedCurrent')
    return
  end function massHaloName

  function massHaloDescription(self)
    !% Return a description of the massHalo property.
    implicit none
    type (varying_string               )                :: massHaloDescription
    class(nodePropertyExtractorMassHalo), intent(inout) :: self
    !$GLC attributes unused :: self

    massHaloDescription=var_str('The current mass of the dark-matter-only halo within a radius enclosing a specified density contrast.')
    return
  end function massHaloDescription

  double precision function massHaloUnitsInSI(self)
    !% Return the units of the massHalo property in the SI system.
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassHalo), intent(inout) :: self
    !$GLC attributes unused :: self

    massHaloUnitsInSI=massSolar
    return
  end function massHaloUnitsInSI
