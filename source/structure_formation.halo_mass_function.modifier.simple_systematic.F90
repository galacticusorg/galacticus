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

!% Contains a module which implements a simple systematic modifier of the dark matter halo mass function.

module Halo_Mass_Function_Modifiers_Simple_Systematic
  !% Implements a simple systematic modifier of the dark matter halo mass function.
  implicit none
  private
  public :: Halo_Mass_Function_Modifier_Simple_Systematic

  ! Record of whether module has been initialized.
  logical          :: moduleInitialized=.false.

  ! Parameters controlling the systematic modifier.
  double precision :: haloMassFunctionSimpleSystematicAlpha,haloMassFunctionSimpleSystematicBeta

contains

  !# <haloMassFunctionModifierMethod>
  !#  <unitName>Halo_Mass_Function_Modifier_Simple_Systematic</unitName>
  !# </haloMassFunctionModifierMethod>
  subroutine Halo_Mass_Function_Modifier_Simple_Systematic(haloTime,haloMass,haloMassFunction)
    !% Modify the halo mass function by applying a simple systematic shift in abundance.
    use Input_Parameters
    implicit none
    double precision, intent(in   )            :: haloTime,haloMass
    double precision, intent(inout)            :: haloMassFunction
    double precision,                parameter :: haloMassZeroPoint=1.0d12

    if (.not.moduleInitialized) then
       !$omp critical (Halo_Mass_Function_Modifier_Simple_Systematic_Initialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>haloMassFunctionSimpleSystematicAlpha</name>
          !@   <defaultValue>$0$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Parameter $\alpha$ appearing in model for simple systematic shift in the halo mass function.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("haloMassFunctionSimpleSystematicAlpha",haloMassFunctionSimpleSystematicAlpha,defaultValue=0.0d0)
          !@ <inputParameter>
          !@   <name>haloMassFunctionSimpleSystematicBeta</name>
          !@   <defaultValue>$0$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Parameter $\beta$ appearing in model for simple systematic shift in the halo mass function.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("haloMassFunctionSimpleSystematicBeta",haloMassFunctionSimpleSystematicBeta,defaultValue=0.0d0)
          ! Record that the module is now initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical (Halo_Mass_Function_Modifier_Simple_Systematic_Initialize)
    end if
    
    ! Compute the systematic shift in the halo mass function.
    haloMassFunction=                               &
         &  haloMassFunction                        &
         & *(                                       &
         &    1.0d0                                 &
         &   +haloMassFunctionSimpleSystematicAlpha &
         &   +haloMassFunctionSimpleSystematicBeta  &
         &   *log10(                                &
         &           haloMass                       &
         &          /haloMassZeroPoint              &
         &         )                                &
         &  )
    return
  end subroutine Halo_Mass_Function_Modifier_Simple_Systematic
  
end module Halo_Mass_Function_Modifiers_Simple_Systematic
