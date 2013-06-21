!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that implements postprocessing of stellar population spectra.

module Stellar_Population_Spectra_Postprocess
  !% Implements postprocessing of stellar population spectra.
  use Abundances_Structure
  use ISO_Varying_String 
  implicit none
  private
  public :: Stellar_Population_Spectrum_Postprocess

  ! Flag to indicate if this module has been initialized.  
  logical                                            :: stellarPopulationSpectraPostprocessInitialized=.false.  
  
  ! Count of the number of methods being applied.                                                                                                           
  integer                                            :: methodCount                                   =0        
  
  ! Name of stellar population postprocessing methods to apply.                                                                                                           
  type   (varying_string), allocatable, dimension(:) :: stellarPopulationSpectraPostprocessMethods              
                                                                                                             
contains

  subroutine Stellar_Population_Spectrum_Postprocess_Initialize
    !% Initialize the stellar population spectra postprocessing module
    use Galacticus_Error
    use Input_Parameters
    use Memory_Management
    !# <include directive="stellarPopulationSpectraPostprocessInitialize" type="moduleUse">
    include 'stellar_populations.spectra.postprocess.initialize.modules.inc'
    !# </include>
    implicit none
    
    ! Initialize if necessary.
    if (.not.stellarPopulationSpectraPostprocessInitialized) then
       !$omp critical(Stellar_Population_Spectrum_Postprocess_Initialization) 
       if (.not.stellarPopulationSpectraPostprocessInitialized) then
          ! Get the stellar population postprocessing methods parameter.
          !@ <inputParameter>
          !@   <name>stellarPopulationSpectraPostprocessMethods</name>
          !@   <defaultValue>Meiksin2006</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of methods to be used for post-processing of stellar population spectra.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          ! Determine how many filters are to be applied.
          methodCount=Get_Input_Parameter_Array_Size('stellarPopulationSpectraPostprocessMethods')
          ! Allocate methods array and read method names.
          if (methodCount > 0) then
             allocate(stellarPopulationSpectraPostprocessMethods(methodCount))
             call Memory_Usage_Record(sizeof(stellarPopulationSpectraPostprocessMethods))
             call Get_Input_Parameter('stellarPopulationSpectraPostprocessMethods',stellarPopulationSpectraPostprocessMethods)
          else
             methodCount=1
             allocate(stellarPopulationSpectraPostprocessMethods(methodCount))
             call Memory_Usage_Record(sizeof(stellarPopulationSpectraPostprocessMethods))
             call Get_Input_Parameter('stellarPopulationSpectraPostprocessMethods',stellarPopulationSpectraPostprocessMethods,defaultValue=['Meiksin2006'])
          end if
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="stellarPopulationSpectraPostprocessInitialize" type="functionCall" functionType="void">
          !#  <functionArgs>stellarPopulationSpectraPostprocessMethods</functionArgs>
          include 'stellar_populations.spectra.postprocess.initialize.inc'
          !# </include>
          stellarPopulationSpectraPostprocessInitialized=.true.
       end if
       !$omp end critical(Stellar_Population_Spectrum_Postprocess_Initialization) 
    end if
    return
  end subroutine Stellar_Population_Spectrum_Postprocess_Initialize

  double precision function Stellar_Population_Spectrum_Postprocess(wavelength,redshift)
    !% Return a multiplicative factor by which a stellar population spectrum should be modified by any postprocessing.
    !# <include directive="stellarPopulationSpectraPostprocess" type="moduleUse">
    include 'stellar_populations.spectra.postprocess.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in   ) :: redshift, wavelength  
    
    ! Initialize the module.                                                      
    call Stellar_Population_Spectrum_Postprocess_Initialize()

    ! Return immediately if no methods were defined.
    if (methodCount == 0) return

    ! Compute the postprocessing factor.
    Stellar_Population_Spectrum_Postprocess=1.0d0
    !# <include directive="stellarPopulationSpectraPostprocess" type="functionCall" functionType="void">
    !#  <functionArgs>wavelength,redshift,Stellar_Population_Spectrum_Postprocess</functionArgs>
    include 'stellar_populations.spectra.postprocess.inc'
    !# </include>
    return
  end function Stellar_Population_Spectrum_Postprocess
  
end module Stellar_Population_Spectra_Postprocess
