!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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


!% Contains a module which implements postprocessing of stellar spectra to keep only recent populations.

module Stellar_Population_Spectra_Postprocessing_Recent
  !% Implements postprocessing of stellar spectra to keep only recent populations.
  use ISO_Varying_String
  public :: Stellar_Population_Spectra_Postprocess_Recent_Init

  ! Record of whether the module is initialized.
  logical          :: moduleInitialized=.false.

  ! Parameters of the model.
  double precision :: recentPopulationsTimeLimit

contains

  !# <stellarPopulationSpectraPostprocessInitialize>
  !#  <unitName>Stellar_Population_Spectra_Postprocess_Recent_Init</unitName>
  !# </stellarPopulationSpectraPostprocessInitialize>
  subroutine Stellar_Population_Spectra_Postprocess_Recent_Init(stellarPopulationSpectraPostprocessMethod,postprocessingFunction)
    !% Initializes the ``recent'' stellar spectrum postprocessing module.
    use Input_Parameters
    implicit none
    type     (varying_string),          intent(in   ) :: stellarPopulationSpectraPostprocessMethod
    procedure(              ), pointer, intent(inout) :: postprocessingFunction

    if (stellarPopulationSpectraPostprocessMethod == 'recent') then
       postprocessingFunction => Stellar_Population_Spectra_Postprocess_Recent
       
       if (.not.moduleInitialized) then
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
          call Get_Input_Parameter('recentPopulationsTimeLimit',recentPopulationsTimeLimit,defaultValue=1.0d-2)
          moduleInitialized=.true.
       end if
    end if
    return
  end subroutine Stellar_Population_Spectra_Postprocess_Recent_Init
  
  subroutine Stellar_Population_Spectra_Postprocess_Recent(wavelength,age,redshift,modifier)
    !% Apply dust attenuation to galaxy spectra
    use Numerical_Constants_Atomic
    implicit none
    double precision, intent(in   ) :: wavelength,age,redshift
    double precision, intent(inout) :: modifier

    if (age > recentPopulationsTimeLimit) modifier=0.0d0
    return
  end subroutine Stellar_Population_Spectra_Postprocess_Recent

end module Stellar_Population_Spectra_Postprocessing_Recent
