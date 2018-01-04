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

  !% An implementation of a spectrum postprocessor that keeps only unescaped populations.

  !# <spectraPostprocessor name="spectraPostprocessorUnescaped">
  !#  <description>Retains only unescaped stellar populations.</description>
  !# </spectraPostprocessor>

  type, extends(spectraPostprocessorClass) :: spectraPostprocessorUnescaped
     !% An unescaped spectrum postprocessor.
     private
     double precision :: timescale
   contains
     procedure :: apply => unescapedApply
  end type spectraPostprocessorUnescaped

  interface spectraPostprocessorUnescaped
     !% Constructors for the unescaped spectrum postprocessor class.
     module procedure unescapedDefaultConstructor
     module procedure unescapedGenericConstructor
  end interface spectraPostprocessorUnescaped

  logical          :: unescapedInitialized=.false.
  double precision :: unescapedTimescale

contains

  function unescapedDefaultConstructor()
    !% Default constructor for the unescaped spectrum postprocessor class.
    use Input_Parameters
    implicit none
    type(spectraPostprocessorUnescaped), target :: unescapedDefaultConstructor
    
    if (.not.unescapedInitialized) then
       ! Get parameters of the model.
       !# <inputParameter>
       !#   <name>stellarPopulationSpectraUnescapedTimescale</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>1.0d-2</defaultValue>
       !#   <description>The timescale for ``escape'' of stellar populations to retain in the ``unescaped'' spectra postprocessing method.</description>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !#   <variable>unescapedTimescale</variable>
       !# </inputParameter>
       unescapedInitialized=.true.
    end if
    unescapedDefaultConstructor=unescapedGenericConstructor(unescapedTimescale)
    return
  end function unescapedDefaultConstructor

  function unescapedGenericConstructor(timescale)
    !% Generic constructor for the unescaped spectrum postprocessor class.
    implicit none
    type            (spectraPostprocessorUnescaped)                :: unescapedGenericConstructor
    double precision                               , intent(in   ) :: timescale

    unescapedGenericConstructor%timescale=timescale
    return
  end function unescapedGenericConstructor

  subroutine unescapedApply(self,wavelength,age,redshift,modifier)
    !% Perform an unescaped postprocessing on a spectrum.
    implicit none
    class           (spectraPostprocessorUnescaped), intent(inout) :: self
    double precision                               , intent(in   ) :: age     , redshift, wavelength
    double precision                               , intent(inout) :: modifier
    !GCC$ attributes unused :: redshift, wavelength
    
    modifier=exp(-age/self%timescale)
    return
  end subroutine unescapedApply
