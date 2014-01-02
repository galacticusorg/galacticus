!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  double precision            :: initialMassGlobal                                      , metallicityGlobal
  !$omp threadprivate(initialMassGlobal,metallicityGlobal)
contains

  !# <stellarFeedbackMethod>
  !#  <unitName>Stellar_Feedback_Standard_Initialize</unitName>
  !# </stellarFeedbackMethod>
  subroutine Stellar_Feedback_Standard_Initialize(stellarFeedbackMethod,Stellar_Feedback_Cumulative_Energy_Input_Get)
    !% Initialize the ``standard'' stellar feedback module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string  ), intent(in   )          :: stellarFeedbackMethod
    procedure(Stellar_Feedback_Cumulative_Energy_Input_Standard), intent(inout), pointer :: Stellar_Feedback_Cumulative_Energy_Input_Get

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
    use Stellar_Astrophysics
    use Supernovae_Type_Ia
    use Supernovae_Population_III
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: age                 , initialMass, metallicity
    double precision                                            :: energySNe           , energyWinds, lifetime
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

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
    real(kind=c_double)        :: Wind_Energy_Integrand
    real(kind=c_double), value :: age
    type(c_ptr        ), value :: parameterPointer

    Wind_Energy_Integrand=0.5d0*Stellar_Winds_Mass_Loss_Rate(initialMassGlobal,age,metallicityGlobal)&
         &*Stellar_Winds_Terminal_Velocity(initialMassGlobal,age,metallicityGlobal)**2
    return
  end function Wind_Energy_Integrand

end module Stellar_Feedback_Standard
