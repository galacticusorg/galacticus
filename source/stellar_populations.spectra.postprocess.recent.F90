!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of a spectrum postprocessor that keeps only recent populations.

  !# <spectraPostprocessor name="spectraPostprocessorRecent">
  !#  <description>Retains only recent stellar populations.</description>
  !# </spectraPostprocessor>

  type, extends(spectraPostprocessorClass) :: spectraPostprocessorRecent
     !% An recent spectrum postprocessor.
     private
     double precision :: timeLimit
   contains
     procedure :: apply => recentApply
  end type spectraPostprocessorRecent

  interface spectraPostprocessorRecent
     !% Constructors for the recent spectrum postprocessor class.
     module procedure recentDefaultConstructor
     module procedure recentGenericConstructor
  end interface spectraPostprocessorRecent

  logical          :: recentInitialized=.false.
  double precision :: recentTimeLimit

contains

  function recentDefaultConstructor()
    !% Default constructor for the recent spectrum postprocessor class.
    use Input_Parameters
    implicit none
    type(spectraPostprocessorRecent), target :: recentDefaultConstructor
    
    if (.not.recentInitialized) then
       ! Get parameters of the model.
       !@ <inputParameter>
       !@   <name>recentPopulationsTimeLimit</name>
       !@   <defaultValue>$10^7$ years</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum age of stellar populations to retain in the ``recent'' spectra postprocessing method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationSpectraRecentTimeLimit',recentTimeLimit,defaultValue=1.0d-2)
       recentInitialized=.true.
    end if
    recentDefaultConstructor=recentGenericConstructor(recentTimeLimit)
    return
  end function recentDefaultConstructor

  function recentGenericConstructor(timeLimit)
    !% Generic constructor for the recent spectrum postprocessor class.
    implicit none
    type            (spectraPostprocessorRecent)                :: recentGenericConstructor
    double precision                            , intent(in   ) :: timeLimit

    recentGenericConstructor%timeLimit=timeLimit
    return
  end function recentGenericConstructor

  subroutine recentApply(self,wavelength,age,redshift,modifier)
    !% Perform an recent postprocessing on a spectrum.
    implicit none
    class           (spectraPostprocessorRecent), intent(inout) :: self
    double precision                            , intent(in   ) :: age     , redshift, wavelength
    double precision                            , intent(inout) :: modifier

    if (age > self%timeLimit) modifier=0.0d0
    return
  end subroutine recentApply
