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

!% Contains a module which implements empirical models of conditional stellar mass functions.

module Conditional_Stellar_Mass_Functions
  !% Implements empirical models of conditional stellar mass functions.
  use ISO_Varying_String
  implicit none
  private
  public :: Cumulative_Conditional_Stellar_Mass_Function,Cumulative_Conditional_Stellar_Mass_Function_Variance
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: conditionalStellarMassFunctionInitialized=.false.

  ! Name of conditional stellar mass function method used.
  type(varying_string) :: conditionalStellarMassFunctionMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Cumulative_Conditional_Stellar_Mass_Function_Template    ), pointer :: Cumulative_Conditional_Stellar_Mass_Function_Get     => null()
  procedure(Cumulative_Conditional_Stellar_Mass_Function_Var_Template), pointer :: Cumulative_Conditional_Stellar_Mass_Function_Var_Get => null()
  abstract interface
     double precision function Cumulative_Conditional_Stellar_Mass_Function_Template(massHalo,massStellar)
       double precision, intent(in) :: massHalo,massStellar
     end function Cumulative_Conditional_Stellar_Mass_Function_Template
  end interface
  abstract interface
     double precision function Cumulative_Conditional_Stellar_Mass_Function_Var_Template(massHalo,massStellarLow,massStellarHigh)
       double precision, intent(in) :: massHalo,massStellarLow,massStellarHigh
     end function Cumulative_Conditional_Stellar_Mass_Function_Var_Template
  end interface

contains

  subroutine Conditional_Stellar_Mass_Functions_Initialize
    !% Initialize the conditional stellar mass function module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="conditionalStellarMassFunctionMethod" type="moduleUse">
    include 'halo_model.conditional_stellar_mass_function.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Conditional_Stellar_Mass_Functions_Initialization) 
    ! Initialize if necessary.
    if (.not.conditionalStellarMassFunctionInitialized) then
       ! Get the conditional stellar mass function method parameter.
       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionMethod</name>
       !@   <defaultValue>Behroozi2010</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for empirical models of the conditional stellar mass function.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionMethod',conditionalStellarMassFunctionMethod,defaultValue='Behroozi2010')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="conditionalStellarMassFunctionMethod" type="functionCall" functionType="void">
       !#  <functionArgs>conditionalStellarMassFunctionMethod,Cumulative_Conditional_Stellar_Mass_Function_Get,Cumulative_Conditional_Stellar_Mass_Function_Var_Get</functionArgs>
       include 'halo_model.conditional_stellar_mass_function.inc'
       !# </include>
       if     (                                                                            &
            &  .not.(                                                                      &
            &             associated(Cumulative_Conditional_Stellar_Mass_Function_Get    ) &
            &        .and.associated(Cumulative_Conditional_Stellar_Mass_Function_Var_Get) &
            &       )                                                                      &
            & ) call Galacticus_Error_Report('Conditional_Stellar_Mass_Functions','method '//char(conditionalStellarMassFunctionMethod)//' is unrecognized')
       conditionalStellarMassFunctionInitialized=.true.
    end if
    !$omp end critical(Conditional_Stellar_Mass_Functions_Initialization) 

    return
  end subroutine Conditional_Stellar_Mass_Functions_Initialize

  double precision function Cumulative_Conditional_Stellar_Mass_Function(massHalo,massStellar)
    !% Returns the cumulative conditional stellar mass function at a stellar mass of {\tt massStellar} in a halo of mass {\tt massHalo}.
    implicit none
    double precision, intent(in) :: massHalo,massStellar
    
    ! Initialize the module.
    call Conditional_Stellar_Mass_Functions_Initialize

    ! Get the mass function using the selected method.
    Cumulative_Conditional_Stellar_Mass_Function=Cumulative_Conditional_Stellar_Mass_Function_Get(massHalo,massStellar)

    return
  end function Cumulative_Conditional_Stellar_Mass_Function
  
  double precision function Cumulative_Conditional_Stellar_Mass_Function_Variance(massHalo,massStellarLow,massStellarHigh)
    !% Returns the cumulative conditional stellar mass function at a stellar mass of {\tt massStellar} in a halo of mass {\tt massHalo}.
    implicit none
    double precision, intent(in) :: massHalo,massStellarLow,massStellarHigh
    
    ! Initialize the module.
    call Conditional_Stellar_Mass_Functions_Initialize

    ! Get the variance using the selected method.
    Cumulative_Conditional_Stellar_Mass_Function_Variance=Cumulative_Conditional_Stellar_Mass_Function_Var_Get(massHalo,massStellarLow,massStellarHigh)

    return
  end function Cumulative_Conditional_Stellar_Mass_Function_Variance
  
end module Conditional_Stellar_Mass_Functions
