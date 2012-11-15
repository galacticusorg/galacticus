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

!% Contains a module which implements halo mass function sampling proportional to halo abundance (i.e. a volume-limited sample of halos).

module Merger_Trees_Mass_Function_Sampling_Halo_MF
  !% Implements halo mass function sampling proportional to halo abundance (i.e. a volume-limited sample of halos).
  private
  public :: Merger_Trees_Mass_Function_Sampling_Halo_MF_Initialize

contains

  !# <haloMassFunctionSamplingMethod>
  !#  <unitName>Merger_Trees_Mass_Function_Sampling_Halo_MF_Initialize</unitName>
  !# </haloMassFunctionSamplingMethod>
  subroutine Merger_Trees_Mass_Function_Sampling_Halo_MF_Initialize(haloMassFunctionSamplingMethod,Merger_Tree_Construct_Mass_Function_Sampling_Get)
    !% Initializes the ``haloMassFunction'' halo mass function sampling method.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: haloMassFunctionSamplingMethod
    procedure(double precision), pointer, intent(inout) :: Merger_Tree_Construct_Mass_Function_Sampling_Get
    
    if (haloMassFunctionSamplingMethod == 'haloMassFunction') Merger_Tree_Construct_Mass_Function_Sampling_Get => Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF
    return
  end subroutine Merger_Trees_Mass_Function_Sampling_Halo_MF_Initialize

  double precision function Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF(mass,time,massMinimum,massMaximum)
    !% Computes the halo mass function sampling rate using a volume-limited sampling.
    use Halo_Mass_Function
    implicit none
    double precision, intent(in) :: mass,time,massMinimum,massMaximum

    Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF=mass*Halo_Mass_Function_Differential(time,mass)
    return
  end function Merger_Tree_Construct_Mass_Function_Sampling_Halo_MF

end module Merger_Trees_Mass_Function_Sampling_Halo_MF
