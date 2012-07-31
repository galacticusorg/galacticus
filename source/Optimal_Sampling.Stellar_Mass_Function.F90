!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu

!% Contains a program which computes the optimal sampling of the halo mass function to minimize errors on the stellar mass
!% function.

program Optimal_Sampling_SMF
  !% Compute the optimal number of trees to run of each mass.
  use, intrinsic :: ISO_C_Binding
  use Memory_Management
  use Numerical_Ranges
  use Input_Parameters
  use ISO_Varying_String
  use FGSL
  use Numerical_Integration
  use Halo_Mass_Function
  use Galacticus_Meta_Compute_Times
  use Cosmology_Functions
  use Merger_Trees_Mass_Function_Sampling
  implicit none
  integer,          parameter                 :: stellarMassPointsPerDecade=100
  integer,          parameter                 ::    haloMassPointsPerDecade=100
  double precision, parameter                 :: stellarMassMinimum=1.0d8 , stellarMassMaximum=1.0d13
  double precision, parameter                 :: haloMassMinimum   =1.0d10, haloMassMaximum   =1.0d15
  double precision, parameter                 :: toleranceAbsolute =0.0d0,  toleranceRelative =1.0d-3
  double precision, allocatable, dimension(:) :: stellarMassTableMass,stellarMassTableMassFunction,haloMassTableMass,haloMassTableXi
  integer                                     :: stellarMassTableCount,haloMassTableCount,iMass,iUnit
  double precision                            :: stellarMass,time,samplingDensity,haloMass
  type(varying_string)                        :: parameterFile
  type(c_ptr)                                 :: parameterPointer
  type(fgsl_function)                         :: integrandFunction
  type(fgsl_integration_workspace)            :: integrationWorkspace
  logical                                     :: integrationReset=.true.

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Optimal_Sampling.Stellar_Mass_Function.size')

  ! Open the parameter file.
  parameterFile='parameters/optimalSamplingStellarMassFunction.xml'
  call Input_Parameters_File_Open(parameterFile)

  ! Compute number of points to tabulate.
  stellarMassTableCount=int(dlog10(stellarMassMaximum/stellarMassMinimum)*dble(stellarMassPointsPerDecade)+1.0d0)
  haloMassTableCount   =int(dlog10(   haloMassMaximum/   haloMassMinimum)*dble(   haloMassPointsPerDecade)+1.0d0)
  
  ! Allocate arrays.
  call Alloc_Array(stellarMassTableMass        ,[stellarMassTableCount])
  call Alloc_Array(stellarMassTableMassFunction,[stellarMassTableCount])
  call Alloc_Array(   haloMassTableMass        ,[   haloMassTableCount])
  call Alloc_Array(   haloMassTableXi          ,[   haloMassTableCount])

  ! Create mass tabulations.
  stellarMassTableMass=Make_Range(stellarMassMinimum,stellarMassMaximum,stellarMassTableCount,rangeType=rangeTypeLogarithmic)
  haloMassTableMass   =Make_Range(   haloMassMinimum,   haloMassMaximum,   haloMassTableCount,rangeType=rangeTypeLogarithmic)

  ! Get the cosmic time at the present day.
  time=Cosmology_Age(1.0d0)

  ! Open a file for the stellar mass function data.
  open(newunit=iUnit,file='stellarMassFunction.data',status='unknown',form='formatted')

  ! Loop over stellar masses.
  do iMass=1,stellarMassTableCount

     ! Set the stellar mass to use.
     stellarMass=stellarMassTableMass(iMass)

     ! Integrate over the mass function weighting by the stellar conditional mass function.
     stellarMassTableMassFunction(iMass)=Integrate(haloMassMinimum,haloMassMaximum,Stellar_Mass_Function_Integrand,parameterPointer&
          &,integrandFunction ,integrationWorkspace,toleranceAbsolute=toleranceAbsolute,toleranceRelative=toleranceRelative,reset&
          &=integrationReset)

     ! Convert to per decade.
     stellarMassTableMassFunction(iMass)=stellarMassTableMassFunction(iMass)*log(10.0d0)

     ! Write out the result.
     write (iUnit,*) iMass,stellarMassTableMass(iMass),stellarMassTableMassFunction(iMass)

  end do
  call Integrate_Done(integrandFunction,integrationWorkspace)
  integrationReset=.true.
  close(iUnit)

   ! Open a file for the tree weighting function data.
   open(newunit=iUnit,file='treeSampling.data',status='unknown',form='formatted')

   ! Loop over halo masses.
   do iMass=1,haloMassTableCount

      ! Set the halo mass to use.
      haloMass=haloMassTableMass(iMass)

      ! Get the sampling density function.
      samplingDensity=Merger_Tree_Construct_Mass_Function_Sampling(haloMass,time,haloMassMinimum,haloMassMaximum)

      ! Write out the result.
      write (iUnit,*) iMass,haloMassTableMass(iMass),samplingDensity,Halo_Mass_Function_Differential(time,haloMass)*haloMass&
           &,Galacticus_Time_Per_Tree(haloMass)
      
   end do
   close(iUnit)

  ! Close the parameter file.
  call Input_Parameters_File_Close

contains

  function Stellar_Mass_Function_Integrand(mass,parameterPointer) bind(c)
    !% The integrand (as a function of halo mass) giving the stellar mass function.
    use, intrinsic :: ISO_C_Binding
    use Halo_Mass_Function
    use Conditional_Stellar_Mass_Functions
    implicit none
    real(c_double)              :: Stellar_Mass_Function_Integrand
    real(c_double),   value     :: mass
    type(c_ptr),      value     :: parameterPointer
    double precision, parameter :: deltaLogMass=0.097d0
    double precision            :: conditionalStellarMassFunction

    conditionalStellarMassFunction=&
         & (Cumulative_Conditional_Stellar_Mass_Function(mass,stellarMass*(10.0d0**(+0.5d0*deltaLogMass))) &
         & -Cumulative_Conditional_Stellar_Mass_Function(mass,stellarMass*(10.0d0**(-0.5d0*deltaLogMass))) &
         & /(10.0d0**(+0.5d0*deltaLogMass)-10.0d0**(-0.5d0*deltaLogMass)))

    Stellar_Mass_Function_Integrand=Halo_Mass_Function_Differential(time,mass)*conditionalStellarMassFunction
    return
  end function Stellar_Mass_Function_Integrand

end program Optimal_Sampling_SMF
