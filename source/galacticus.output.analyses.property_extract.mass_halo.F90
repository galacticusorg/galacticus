!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  use Virial_Density_Contrast

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorMassHalo">
  !#  <description>A halo mass output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorMassHalo
     !% A halo mass property extractor output analysis class. The property extracted is the ''\gls{dmou}'' mass of the halo within
     !% a radius enclosing a density contrast as defined by the supplied {\normalfont \ttfamily virialDensityContrast} class
     !% object. Note that the density contrast is defined here at the time at which the halo presently exists, \emph{not} at the
     !% time at which is was last isolated (as is used for standard definition of halo mass).
     private
     class(virialDensityContrastClass), pointer :: virialDensityContrast_
   contains
     procedure :: extract  => massHaloExtract
     procedure :: type     => massHaloType
  end type outputAnalysisPropertyExtractorMassHalo

  interface outputAnalysisPropertyExtractorMassHalo
     !% Constructors for the ``massHalo'' output analysis class.
     module procedure massHaloConstructorParameters
     module procedure massHaloConstructorInternal
  end interface outputAnalysisPropertyExtractorMassHalo

contains

  function massHaloConstructorParameters(parameters) result(self)
    !% Constructor for the ``massHalo'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyExtractorMassHalo)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(virialDensityContrastClass             ), pointer       :: virialDensityContrast_
    !GCC$ attributes unused :: parameters

    ! Get the default virial density contrast object.
    virialDensityContrast_ => virialDensityContrast()
    ! Build the object.
    self=outputAnalysisPropertyExtractorMassHalo(virialDensityContrast_)
    return
  end function massHaloConstructorParameters

  function massHaloConstructorInternal(virialDensityContrast_) result(self)
    !% Internal constructor for the ``massHalo'' output analysis property extractor class.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyExtractorMassHalo)         :: self
    class(virialDensityContrastClass             ), target :: virialDensityContrast_
    !# <constructorAssign variables="*virialDensityContrast_"/>
    
    return
  end function massHaloConstructorInternal

  double precision function massHaloExtract(self,node)
    !% Implement a massHalo output analysis.
    use Dark_Matter_Profile_Mass_Definitions
    implicit none
    class(outputAnalysisPropertyExtractorMassHalo), intent(inout) :: self
    type (treeNode                               ), intent(inout) :: node
    class(nodeComponentBasic                     ), pointer       :: basic
    
    basic           => node%basic()
    massHaloExtract =  Dark_Matter_Profile_Mass_Definition(node,self%virialDensityContrast_%densityContrast(basic%mass(),basic%time()))    
    return
  end function massHaloExtract

  integer function massHaloType(self)
    !% Return the type of the halo mass property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassHalo), intent(inout) :: self
    !GCC$ attributes unused :: self

    massHaloType=outputAnalysisPropertyTypeLinear
    return
  end function massHaloType
