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

!% Contains a module which implements halo mass function sampling using a power-law in halo mass.

module Merger_Trees_Mass_Function_Sampling_Power_Law
  !% Implements halo mass function sampling using a power-law in halo mass.
  private
  public :: Merger_Trees_Mass_Function_Sampling_Power_Law_Initialize

  ! The exponent in the halo mass function sampling power law.
  double precision :: mergerTreeBuildTreesHaloMassExponent

contains

  !# <haloMassFunctionSamplingMethod>
  !#  <unitName>Merger_Trees_Mass_Function_Sampling_Power_Law_Initialize</unitName>
  !# </haloMassFunctionSamplingMethod>
  subroutine Merger_Trees_Mass_Function_Sampling_Power_Law_Initialize(haloMassFunctionSamplingMethod,Merger_Tree_Construct_Mass_Function_Sampling_Get)
    !% Initializes the ``powerLaw'' halo mass function sampling method.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: haloMassFunctionSamplingMethod
    procedure(double precision), pointer, intent(inout) :: Merger_Tree_Construct_Mass_Function_Sampling_Get
    
    if (haloMassFunctionSamplingMethod == 'powerLaw') then
       Merger_Tree_Construct_Mass_Function_Sampling_Get => Merger_Tree_Construct_Mass_Function_Sampling_Power_Law

       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesHaloMassExponent</name>
       !@   <defaultValue>1</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Halo masses will be (pseudo-)uniformly distributed in $[\log(M)]^{1/(1+\alpha)}$ where $\alpha=${\tt mergerTreeBuildTreesHaloMassExponent}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesHaloMassExponent',mergerTreeBuildTreesHaloMassExponent,defaultValue=1.0d0)

    end if
    return
  end subroutine Merger_Trees_Mass_Function_Sampling_Power_Law_Initialize

  double precision function Merger_Tree_Construct_Mass_Function_Sampling_Power_Law(mass,time,massMinimum,massMaximum)
    !% Computes the halo mass function sampling rate using a power-law distribution.
    implicit none
    double precision, intent(in) :: mass,time,massMinimum,massMaximum
    
    ! Sampling rate is simply a power-law in the logarithm of halo mass.
    if (mass <= massMinimum .or. mass > massMaximum) then
       Merger_Tree_Construct_Mass_Function_Sampling_Power_Law=0.0d0
    else
       Merger_Tree_Construct_Mass_Function_Sampling_Power_Law=1.0d0/log10(mass/massMinimum)**(mergerTreeBuildTreesHaloMassExponent/(1.0d0+mergerTreeBuildTreesHaloMassExponent))
    end if
    return
  end function Merger_Tree_Construct_Mass_Function_Sampling_Power_Law

end module Merger_Trees_Mass_Function_Sampling_Power_Law
