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

!% Contains a module which implements a radiusHalfMass property extractor class.

  !# <nodePropertyExtractor name="nodePropertyExtractorRadiusHalfMass">
  !#  <description>
  !#   A node property extractor which extracts the half-mass radius of the galaxy. The half-mass radius is output as {\normalfont
  !#   \ttfamily [halfMassRadius]} (in Mpc).
  !#  </description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusHalfMass
     !% A half-mass radius property extractor class.
     private
   contains
     procedure :: extract     => radiusHalfMassExtract
     procedure :: name        => radiusHalfMassName
     procedure :: description => radiusHalfMassDescription
     procedure :: unitsInSI   => radiusHalfMassUnitsInSI
     procedure :: type        => radiusHalfMassType
  end type nodePropertyExtractorRadiusHalfMass

  interface nodePropertyExtractorRadiusHalfMass
     !% Constructors for the ``radiusHalfMass'' output analysis class.
     module procedure radiusHalfMassConstructorParameters
  end interface nodePropertyExtractorRadiusHalfMass

contains

  function radiusHalfMassConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily radiusHalfMass} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiusHalfMass)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=nodePropertyExtractorRadiusHalfMass()
    return
  end function radiusHalfMassConstructorParameters

  double precision function radiusHalfMassExtract(self,node,instance)
    !% Implement a last isolated redshift output analysis.
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Radius_Enclosing_Mass
    use :: Galactic_Structure_Options        , only : massTypeStellar
    implicit none
    class(nodePropertyExtractorRadiusHalfMass), intent(inout)           :: self
    type (treeNode                           ), intent(inout), target   :: node
    type (multiCounter                       ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance

    radiusHalfMassExtract=Galactic_Structure_Radius_Enclosing_Mass(node,fractionalMass=0.5d0,massType=massTypeStellar)
    return
  end function radiusHalfMassExtract

  function radiusHalfMassName(self)
    !% Return the name of the last isolated redshift property.
    implicit none
    type (varying_string                     )                :: radiusHalfMassName
    class(nodePropertyExtractorRadiusHalfMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassName=var_str('radiusHalfMassStellar')
    return
  end function radiusHalfMassName

  function radiusHalfMassDescription(self)
    !% Return a description of the radiusHalfMass property.
    implicit none
    type (varying_string                     )                :: radiusHalfMassDescription
    class(nodePropertyExtractorRadiusHalfMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassDescription=var_str('Radius enclosing half the galaxy stellar mass [Mpc].')
    return
  end function radiusHalfMassDescription

  double precision function radiusHalfMassUnitsInSI(self)
    !% Return the units of the last isolated redshift property in the SI system.
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorRadiusHalfMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassUnitsInSI=massSolar
    return
  end function radiusHalfMassUnitsInSI

  integer function radiusHalfMassType(self)
    !% Return the type of the last isolated redshift property.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorRadiusHalfMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusHalfMassType=outputAnalysisPropertyTypeLinear
    return
  end function radiusHalfMassType
