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

!% Contains a module which implements the Chabrier stellar initial mass function \citep{chabrier_galactic_2001}.

module Star_Formation_IMF_Chabrier
  !% Implements the Chabrier stellar initial mass function.
  implicit none
  private
  public :: Star_Formation_IMF_Register_Chabrier, Star_Formation_IMF_Register_Name_Chabrier,&
       & Star_Formation_IMF_Recycled_Instantaneous_Chabrier, Star_Formation_IMF_Yield_Instantaneous_Chabrier,&
       & Star_Formation_IMF_Tabulate_Chabrier, Star_Formation_IMF_Minimum_Mass_Chabrier, Star_Formation_IMF_Maximum_Mass_Chabrier&
       &, Star_Formation_IMF_Phi_Chabrier

  ! Index assigned to this IMF.
  integer                     :: imfIndex                        =-1

  ! Flag indicating if the module has been initialized.
  logical                     :: imfChabrierInitialized          =.false.

  ! Parameters of the IMF.
  double precision            :: imfChabrierRecycledInstantaneous               , imfChabrierYieldInstantaneous

  ! Fixed parameters of the IMF.
  double precision, parameter :: chabrierMassLower               =0.1d0         , chabrierMassTransition       =1.0d0          , &
       &                         chabrierMassUpper               =125.0d0
  double precision, parameter :: chabrierExponent                =-2.3d0        , chabrierMass                 =0.08d0         , &
       &                         chabrierSigma                   =0.69d0
  double precision, parameter :: normalizationExponential        =0.2353362401d0, normalizationLogNormal       =0.8326617874d0

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_Chabrier</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_Chabrier(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex=imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_Chabrier

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_Chabrier</unitName>
  !#  <name>Chabrier</name>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_Chabrier(imfNames,imfDescriptors)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(inout) :: imfDescriptors(:), imfNames(:)

    imfNames      (imfIndex)="Chabrier"
    imfDescriptors(imfIndex)="Chabrier"
    return
  end subroutine Star_Formation_IMF_Register_Name_Chabrier

  subroutine Star_Formation_IMF_Initialize_Chabrier
    !% Initialize the Chabrier IMF module.
    use Input_Parameters
    implicit none

    if (.not.imfChabrierInitialized) then
       !$omp critical (IMF_Chabrier_Initialize)
       if (.not.imfChabrierInitialized) then
          !@ <inputParameter>
          !@   <name>imfChabrierRecycledInstantaneous</name>
          !@   <defaultValue>0.46 (internally computed)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The recycled fraction for the Chabrier \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierRecycledInstantaneous',imfChabrierRecycledInstantaneous,defaultValue=0.46d0)
          !@ <inputParameter>
          !@   <name>imfChabrierYieldInstantaneous</name>
          !@   <defaultValue>0.035 (internally computed)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The yield for the Chabrier \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierYieldInstantaneous'   ,imfChabrierYieldInstantaneous   ,defaultValue=0.035d0)
          imfChabrierInitialized=.true.
       end if
       !$omp end critical (IMF_Chabrier_Initialize)
    end if
    return
  end subroutine Star_Formation_IMF_Initialize_Chabrier

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_Chabrier</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_Chabrier(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       minimumMass=chabrierMassLower
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_Chabrier

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_Chabrier</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_Chabrier(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       maximumMass=chabrierMassUpper
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_Chabrier

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_Chabrier</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_Chabrier(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(in   ) :: initialMass
    double precision, intent(  out) :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       imfPhi=Chabrier_Phi(initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_Chabrier

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_Chabrier</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_Chabrier(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       recycledFraction=imfChabrierRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_Chabrier

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_Chabrier</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_Chabrier(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       yield=imfChabrierYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_Chabrier

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_Chabrier</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_Chabrier(imfSelected,imfMatched,imf)
    !% Register the name of this IMF.
    use Tables
    implicit none
    integer                      , intent(in   ) :: imfSelected
    logical                      , intent(inout) :: imfMatched
    class  (table1D), allocatable, intent(inout) :: imf
    integer         , parameter                  :: nPoints    =100

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier()
       allocate(table1DLogarithmicLinear :: imf)
       select type (imf)
       type is (table1DLogarithmicLinear)
          call imf%create  (                   &
               &            chabrierMassLower, &
               &            chabrierMassUpper, &
               &            nPoints            &
               &           )
          call imf%populate(Chabrier_Phi(imf%xs()))
          imfMatched=.true.
       end select
    end if
    return
  end subroutine Star_Formation_IMF_Tabulate_Chabrier

  elemental double precision function Chabrier_Phi(mass)
    !% Evaluates the Chabrier initial mass function.
    implicit none
    double precision, intent(in   ) :: mass

    if (mass >= chabrierMassLower .and. mass < chabrierMassTransition) then
       Chabrier_Phi=normalizationLogNormal*exp(-0.5d0*(log10(mass/chabrierMass)/chabrierSigma)**2)/mass
    else if (mass >= chabrierMassTransition .and. mass < chabrierMassUpper) then
       Chabrier_Phi=normalizationExponential*(mass**chabrierExponent)
    else
       Chabrier_Phi=0.0d0
    end if
    return
  end function Chabrier_Phi

end module Star_Formation_IMF_Chabrier
