!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements a simple calculation of energy feedback from stellar populations.

module Stellar_Feedback_Standard
  !% Implements a simple calculation of energy feedback from stellar populations.
  use Numerical_Constants_Astronomical
  implicit none
  private
  public :: Stellar_Feedback_Standard_Initialize

  ! Parameters controlling the module.
  double precision            :: initialMassForSupernovaeTypeII
  double precision            :: supernovaEnergy
  double precision, parameter :: populationIIIMaximumMetallicity=1.0d-4*metallicitySolar

  ! Global variables used in integrands.
  double precision :: initialMassGlobal,metallicityGlobal
  !$omp threadprivate(initialMassGlobal,metallicityGlobal)

contains

  !# <stellarFeedbackMethod>
  !#  <unitName>Stellar_Feedback_Standard_Initialize</unitName>
  !# </stellarFeedbackMethod>
  subroutine Stellar_Feedback_Standard_Initialize(stellarFeedbackMethod,Stellar_Feedback_Cumulative_Energy_Input_Get)
    !% Initialize the ``standard'' stellar feedback module.
    use Numerical_Constants_Units
    use Numerical_Constants_Prefixes
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: stellarFeedbackMethod
    procedure(double precision), pointer, intent(inout) :: Stellar_Feedback_Cumulative_Energy_Input_Get

    if (stellarFeedbackMethod == 'standard') then
       ! Set procedure pointers.
       Stellar_Feedback_Cumulative_Energy_Input_Get => Stellar_Feedback_Cumulative_Energy_Input_Standard

       ! Read in parameters required by this module.
       !@ <inputParameter>
       !@   <name>initialMassForSupernovaeTypeII</name>
       !@   <defaultValue>8</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum mass that a star must have in order that is result in a Type II supernova.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('initialMassForSupernovaeTypeII',initialMassForSupernovaeTypeII,defaultValue=8.0d0)
       !@ <inputParameter>
       !@   <name>supernovaEnergy</name>
       !@   <defaultValue>$10^{51}$ ergs</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The energy produced by a supernova (in ergs).
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('supernovaEnergy',supernovaEnergy,defaultValue=1.0d51)
       ! Convert energy to MSolar (km/s)^2.
       supernovaEnergy=supernovaEnergy*ergs/massSolar/kilo**2
    end if
    return
  end subroutine Stellar_Feedback_Standard_Initialize

  double precision function Stellar_Feedback_Cumulative_Energy_Input_Standard(initialMass,age,metallicity)
    !% Compute the cumulative energy input from a star of given {\tt initialMass}, {\tt age} and {\tt metallicity}.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Stellar_Astrophysics
    use Supernovae_Type_Ia
    use Supernovae_Population_III
    use Numerical_Integration
    implicit none
    double precision,                intent(in) :: initialMass,age,metallicity
    double precision                            :: lifetime,energyWinds,energySNe
    type(c_ptr)                                 :: parameterPointer
    type(fgsl_function)                         :: integrandFunction
    type(fgsl_integration_workspace)            :: integrationWorkspace    

    ! Begin with zero energy input.
    Stellar_Feedback_Cumulative_Energy_Input_Standard=0.0d0

    ! Check if the star is sufficiently massive to result in a Type II supernova.
    if (initialMass >= initialMassForSupernovaeTypeII ) then
       ! Get the lifetime of the star.
       lifetime=Star_Lifetime(initialMass,metallicity)
       ! If lifetime is exceeded, assume a SNe has occurred.
       if (age >= lifetime) then
          energySNe=0.0d0
          ! Check for pair instability supernovae.
          if (metallicity <= populationIIIMaximumMetallicity) energySNe=energySNe+SNePopIII_Cumulative_Energy(initialMass,age&
               &,metallicity)
          ! Population II star - normal supernova.
          energySNe=energySNe+supernovaEnergy
          ! Add the supernova energy.
          Stellar_Feedback_Cumulative_Energy_Input_Standard=Stellar_Feedback_Cumulative_Energy_Input_Standard+energySNe
       end if
    end if

    ! Add in contribution from Type Ia supernovae.
    Stellar_Feedback_Cumulative_Energy_Input_Standard=Stellar_Feedback_Cumulative_Energy_Input_Standard&
         &+SNeIa_Cumulative_Number(initialMass,age,metallicity)*supernovaEnergy

    ! Add in the contribution from stellar winds.
    initialMassGlobal=initialMass
    metallicityGlobal=metallicity
    energyWinds=Integrate(0.0d0,age,Wind_Energy_Integrand,parameterPointer,integrandFunction,integrationWorkspace&
         &,toleranceAbsolute=1.0d-3*Stellar_Feedback_Cumulative_Energy_Input_Standard,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    Stellar_Feedback_Cumulative_Energy_Input_Standard=Stellar_Feedback_Cumulative_Energy_Input_Standard+energyWinds

    return
  end function Stellar_Feedback_Cumulative_Energy_Input_Standard

  function Wind_Energy_Integrand(age,parameterPointer) bind(c)
    !% Integrand used in evaluating cumulative energy input from winds.
    use, intrinsic :: ISO_C_Binding
    use Stellar_Astrophysics_Winds
    implicit none
    real(c_double)        :: Wind_Energy_Integrand
    real(c_double), value :: age
    type(c_ptr),    value :: parameterPointer

    Wind_Energy_Integrand=0.5d0*Stellar_Winds_Mass_Loss_Rate(initialMassGlobal,age,metallicityGlobal)&
         &*Stellar_Winds_Terminal_Velocity(initialMassGlobal,age,metallicityGlobal)**2
    return
  end function Wind_Energy_Integrand
  
end module Stellar_Feedback_Standard
