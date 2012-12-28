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

!% Contains a module which implements methods for sampling the halo mass function when constructing merger trees.

module Merger_Trees_Mass_Function_Sampling
  !% Implements methods for sampling the halo mass function when constructing merger trees.
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Construct_Mass_Function_Sampling
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: haloMassFunctionSamplingInitialized=.false.

  ! Name of conditional stellar mass function method used.
  type(varying_string) :: haloMassFunctionSamplingMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Merger_Tree_Construct_Mass_Function_Sampling_Template), pointer :: Merger_Tree_Construct_Mass_Function_Sampling_Get => null()
  abstract interface
     double precision function Merger_Tree_Construct_Mass_Function_Sampling_Template(mass,time,massMinimum,massMaximum)
       double precision, intent(in) :: mass,time,massMinimum,massMaximum
     end function Merger_Tree_Construct_Mass_Function_Sampling_Template
  end interface

contains

  subroutine Merger_Trees_Mass_Function_Sampling_Initialize
    !% Initialize the halo mass function sampling module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="haloMassFunctionSamplingMethod" type="moduleUse">
    include 'merger_trees.construct.mass_function_sampling.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Merger_Trees_Mass_Function_Sampling_Initialization) 
    ! Initialize if necessary.
    if (.not.haloMassFunctionSamplingInitialized) then
       ! Get the conditional stellar mass function method parameter.
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingMethod</name>
       !@   <defaultValue>powerLaw</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for sampling the halo mass function when constructing merger trees.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingMethod',haloMassFunctionSamplingMethod,defaultValue='powerLaw')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="haloMassFunctionSamplingMethod" type="functionCall" functionType="void">
       !#  <functionArgs>haloMassFunctionSamplingMethod,Merger_Tree_Construct_Mass_Function_Sampling_Get</functionArgs>
       include 'merger_trees.construct.mass_function_sampling.inc'
       !# </include>
       if     (                                                                   &
            &  .not.(                                                             &
            &        associated(Merger_Tree_Construct_Mass_Function_Sampling_Get) &
            &       )                                                             &
            & ) call Galacticus_Error_Report('Merger_Trees_Mass_Function_Sampling','method '//char(haloMassFunctionSamplingMethod)//' is unrecognized')
       haloMassFunctionSamplingInitialized=.true.
    end if
    !$omp end critical(Merger_Trees_Mass_Function_Sampling_Initialization) 

    return
  end subroutine Merger_Trees_Mass_Function_Sampling_Initialize

  double precision function Merger_Tree_Construct_Mass_Function_Sampling(mass,time,massMinimum,massMaximum)
    !% Returns the sampling rate for merger trees of the given {\tt mass}, per decade of halo mass.
    implicit none
    double precision, intent(in) :: mass,time,massMinimum,massMaximum
    
    ! Initialize the module.
    call Merger_Trees_Mass_Function_Sampling_Initialize

    ! Get the sampling rate using the selected method.
    Merger_Tree_Construct_Mass_Function_Sampling=Merger_Tree_Construct_Mass_Function_Sampling_Get(mass,time,massMinimum,massMaximum)

    return
  end function Merger_Tree_Construct_Mass_Function_Sampling

end module Merger_Trees_Mass_Function_Sampling
