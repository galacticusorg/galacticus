!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements a concentration output analysis property extractor class.

  use Virial_Density_Contrast

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorConcentration" defaultThreadPrivate="yes">
  !#  <description>A concentration output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorConcentration
     !% A concentration property extractor output analysis class. The property extracted is the ''\gls{dmou}'' concentration of
     !% the halo with a radius defined by a density contrast as given by the supplied {\normalfont \ttfamily
     !% virialDensityContrast} class object. Note that the density contrast is defined here at the time at which the halo
     !% presently exists, \emph{not} at the time at which is was last isolated (as is used for standard definition of
     !% concentration).
     private
     class(virialDensityContrastClass), pointer :: virialDensityContrast_
   contains
     final     ::             concentrationDestructor
     procedure :: extract  => concentrationExtract
     procedure :: type     => concentrationType
  end type outputAnalysisPropertyExtractorConcentration

  interface outputAnalysisPropertyExtractorConcentration
     !% Constructors for the ``concentration'' output analysis class.
     module procedure concentrationConstructorParameters
     module procedure concentrationConstructorInternal
  end interface outputAnalysisPropertyExtractorConcentration

contains

  function concentrationConstructorParameters(parameters) result(self)
    !% Constructor for the ``concentration'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyExtractorConcentration)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(virialDensityContrastClass                  ), pointer       :: virialDensityContrast_

    !# <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    self=outputAnalysisPropertyExtractorConcentration(virialDensityContrast_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function concentrationConstructorParameters

  function concentrationConstructorInternal(virialDensityContrast_) result(self)
    !% Internal constructor for the ``concentration'' output analysis property extractor class.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyExtractorConcentration)         :: self
    class(virialDensityContrastClass                  ), target :: virialDensityContrast_
    !# <constructorAssign variables="*virialDensityContrast_"/>
    
    return
  end function concentrationConstructorInternal

  subroutine concentrationDestructor(self)
    !% Destructor for the {\normalfont \ttfamily concentration} output analysis property extractor class.
    implicit none
    type(outputAnalysisPropertyExtractorConcentration), intent(inout) :: self

    !# <objectDestructor name="self%virialDensityContrast_"/>
    return
  end subroutine concentrationDestructor
  
  double precision function concentrationExtract(self,node)
    !% Implement a concentration output analysis.
    use Dark_Matter_Profile_Mass_Definitions
    use Galacticus_Nodes                    , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (outputAnalysisPropertyExtractorConcentration), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    class           (nodeComponentBasic                          ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile              ), pointer       :: darkMatterProfile
    double precision                                                              :: massHalo         , radiusHalo

    basic                =>  node%basic            ()
    darkMatterProfile    =>  node%darkMatterProfile()
    massHalo             =   Dark_Matter_Profile_Mass_Definition(                                                                                     &
         &                                                              node                                                                        , &
         &                                                              self      %virialDensityContrast_%densityContrast(basic%mass(),basic%time()), &
         &                                                       radius=radiusHalo                                                                    &
         &                                                      )
    concentrationExtract =  +                  radiusHalo   &
         &                  /darkMatterProfile%scale     ()
    return
  end function concentrationExtract

  integer function concentrationType(self)
    !% Return the type of the concentration property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorConcentration), intent(inout) :: self
    !GCC$ attributes unused :: self

    concentrationType=outputAnalysisPropertyTypeLinear
    return
  end function concentrationType
