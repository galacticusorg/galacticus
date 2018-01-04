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

!% Contains a module which implements an ISM metallicity output analysis property extractor class.

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorMetallicityISM">
  !#  <description>An ISM metallicity output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorMetallicityISM
     !% An ISM metallicity output analysis property extractor class.
     private
   contains
     procedure :: extract => metallicityISMExtract
     procedure :: type    => metallicityISMType
  end type outputAnalysisPropertyExtractorMetallicityISM

  interface outputAnalysisPropertyExtractorMetallicityISM
     !% Constructors for the ``metallicityISM'' output analysis class.
     module procedure metallicityISMConstructorParameters
  end interface outputAnalysisPropertyExtractorMetallicityISM

contains

  function metallicityISMConstructorParameters(parameters) result(self)
    !% Constructor for the ``metallicityISM'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(outputAnalysisPropertyExtractorMetallicityISM)                :: self
    type(inputParameters                              ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=outputAnalysisPropertyExtractorMetallicityISM()
    return
  end function metallicityISMConstructorParameters

  double precision function metallicityISMExtract(self,node)
    !% Implement a metallicityISM output analysis.
    use Abundances_Structure
    implicit none
    class           (outputAnalysisPropertyExtractorMetallicityISM), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    class           (nodeComponentDisk                            ), pointer       :: disk
    class           (nodeComponentSpheroid                        ), pointer       :: spheroid
    type            (abundances                                   )                :: abundancesDisk, abundancesSpheroid
    double precision                                                               :: massGas       , massMetals
    !GCC$ attributes unused :: self

    disk               =>  node          %disk         ()
    spheroid           =>  node          %spheroid     ()
    abundancesDisk     =   disk          %abundancesGas()
    abundancesSpheroid =                                  spheroid          %abundancesGas()
    massGas            =  +disk          %massGas      ()+spheroid          %massGas      ()
    massMetals         =  +abundancesDisk%metallicity  ()+abundancesSpheroid%metallicity  ()
    if (massGas > 0.0d0) then
       metallicityISMExtract=+massMetals &
            &                /massGas
    else
       metallicityISMExtract=0.0d0
    end if
    return
  end function metallicityISMExtract

  integer function metallicityISMType(self)
    !% Return the type of the stellar mass property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorMetallicityISM), intent(inout) :: self
    !GCC$ attributes unused :: self

    metallicityISMType=outputAnalysisPropertyTypeLinear
    return
  end function metallicityISMType
