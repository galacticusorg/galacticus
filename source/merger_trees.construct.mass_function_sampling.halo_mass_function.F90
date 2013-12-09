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

!% Contains a module which implements halo mass function sampling proportional to halo abundance (i.e. a volume-limited sample of halos).

module Merger_Trees_Mass_Function_Sampling_Halo_MF
  !% Implements halo mass function sampling proportional to halo abundance (i.e. a volume-limited sample of halos).
  private
  public :: Merger_Trees_Mass_Function_Sampling_Halo_MF_Initialize

  ! Limits on abundance.
  double precision :: haloMassFunctionSamplingAbundanceMinimum, haloMassFunctionSamplingAbundanceMaximum

contains

  !# <haloMassFunctionSamplingMethod>
  !#  <unitName>Merger_Trees_Mass_Function_Sampling_Halo_MF_Initialize</unitName>
  !# </haloMassFunctionSamplingMethod>
  subroutine Merger_Trees_Mass_Function_Sampling_Halo_MF_Initialize(haloMassFunctionSamplingMethod,Merger_Tree_Construct_Mass_Function_Sampling_Get)
    !% Initializes the ``haloMassFunction'' halo mass function sampling method.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string  ), intent(in   )          :: haloMassFunctionSamplingMethod
    procedure(Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF), intent(inout), pointer :: Merger_Tree_Construct_Mass_Function_Sampling_Get

    if (haloMassFunctionSamplingMethod == 'haloMassFunction') then
       Merger_Tree_Construct_Mass_Function_Sampling_Get => Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingAbundanceMinimum</name>
       !@   <defaultValue>-1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The abundance (in units of Mpc$^{-3}$) below which to truncate the halo mass function when sampling halo masses for tree construction. A negative value indicates no truncation.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingAbundanceMinimum',haloMassFunctionSamplingAbundanceMinimum,defaultValue=-1.0d0)
       !@ <inputParameter>
       !@   <name>haloMassFunctionSamplingAbundanceMaximum</name>
       !@   <defaultValue>-1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The abundance (in units of Mpc$^{-3}$) above which to truncate the halo mass function when sampling halo masses for tree construction. A negative value indicates no truncation.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloMassFunctionSamplingAbundanceMaximum',haloMassFunctionSamplingAbundanceMaximum,defaultValue=-1.0d0)
    end if
    return
  end subroutine Merger_Trees_Mass_Function_Sampling_Halo_MF_Initialize

  double precision function Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF(mass,time,massMinimum,massMaximum)
    !% Computes the halo mass function sampling rate using a volume-limited sampling.
    use Halo_Mass_Function
    implicit none
    double precision, intent(in   ) :: mass, massMaximum, massMinimum, time

    Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF=mass*Halo_Mass_Function_Differential(time,mass)
    if (haloMassFunctionSamplingAbundanceMinimum > 0.0d0)             &
         & Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF       &
         & =max(                                                      &
         &      Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF, &
         &      haloMassFunctionSamplingAbundanceMinimum              &
         &     )
    if (haloMassFunctionSamplingAbundanceMaximum > 0.0d0)             &
         & Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF       &
         & =min(                                                      &
         &      Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF, &
         &      haloMassFunctionSamplingAbundanceMaximum              &
         &     )
   return
  end function Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF

end module Merger_Trees_Mass_Function_Sampling_Halo_MF
