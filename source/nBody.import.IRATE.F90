!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements an N-body data importer for IRATE files.

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: Cosmology_Functions , only : cosmologyFunctionsClass

  !# <nbodyImporter name="nbodyImporterIRATE">
  !#  <description>An importer for IRATE files.</description>
  !# </nbodyImporter>
  type, extends(nbodyImporterClass) :: nbodyImporterIRATE
     !% An importer for IRATE files.
     private
     class  (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class  (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     type   (varying_string          )          :: fileName
     integer                                    :: snapshot
   contains
     final     ::           irateDestructor
     procedure :: import => irateImport
  end type nbodyImporterIRATE

  interface nbodyImporterIRATE
     !% Constructors for the {\normalfont \ttfamily irate} N-body importer class.
     module procedure irateConstructorParameters
     module procedure irateConstructorInternal
  end interface nbodyImporterIRATE

contains

  function irateConstructorParameters(parameters) result (self)
    !% Constructor for the {\normalfont \ttfamily irate} N-body importer class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyImporterIRATE      )                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    class  (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class  (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    type   (varying_string          )                :: fileName
    integer                                          :: snapshot

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of the file to read.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>snapshot</name>
    !#   <source>parameters</source>
    !#   <description>The snapshot number to read.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    self=nbodyImporterIRATE(fileName,snapshot,cosmologyParameters_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    return
  end function irateConstructorParameters

  function irateConstructorInternal(fileName,snapshot,cosmologyParameters_,cosmologyFunctions_) result (self)
    !% Internal constructor for the {\normalfont \ttfamily irate} N-body importer class.
    implicit none
    type   (nbodyImporterIRATE      )                        :: self
    type   (varying_string          ), intent(in   )         :: fileName
    integer                          , intent(in   )         :: snapshot
    class  (cosmologyParametersClass), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="fileName, snapshot, *cosmologyParameters_, *cosmologyFunctions_"/>

    return
  end function irateConstructorInternal

  subroutine irateDestructor(self)
    !% Destructor for {\normalfont \ttfamily irate} importer class.
    implicit none
    type(nbodyImporterIRATE), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%cosmologyFunctions_" />
    return
  end subroutine irateDestructor

  function irateImport(self) result(simulation)
    !% Import data from a IRATE file.
    use :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, verbosityStandard
    use :: Hashes            , only : rank1IntegerSizeTHash    , rank1DoubleHash
    use :: IO_IRATE          , only : irate
    implicit none
    type (nBodyData         )                :: simulation
    class(nbodyImporterIRATE), intent(inout) :: self
    type (irate             )                :: irate_

    call Galacticus_Display_Indent('import simulation from IRATE file',verbosityStandard)
    irate_=irate(char(self%fileName),self%cosmologyParameters_,self%cosmologyFunctions_)
    call irate_%readHalos(self%snapshot,center=simulation%position,velocity=simulation%velocity,IDs=simulation%particleIDs)
    simulation%propertiesInteger=rank1IntegerSizeTHash()
    simulation%propertiesReal   =rank1DoubleHash      ()
    call Galacticus_Display_Unindent('done',verbosityStandard)
    return
  end function irateImport
