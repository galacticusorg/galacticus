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

!% Contains a module which implements the Salpeter stellar initial mass function \citep{salpeter_luminosity_1955}.

module Star_Formation_IMF_Salpeter
  !% Implements the Salpeter stellar initial mass function.
  implicit none
  private
  public :: Star_Formation_IMF_Register_Salpeter, Star_Formation_IMF_Register_Name_Salpeter,&
       & Star_Formation_IMF_Recycled_Instantaneous_Salpeter, Star_Formation_IMF_Yield_Instantaneous_Salpeter,&
       & Star_Formation_IMF_Tabulate_Salpeter, Star_Formation_IMF_Minimum_Mass_Salpeter, Star_Formation_IMF_Maximum_Mass_Salpeter&
       &, Star_Formation_IMF_Phi_Salpeter

  ! Index assigned to this IMF.
  integer :: imfIndex=-1

  ! Flag indicating if the module has been initialized.
  logical :: imfSalpeterInitialized=.false.

  ! Parameters of the IMF.
  double precision :: imfSalpeterRecycledInstantaneous, imfSalpeterYieldInstantaneous

  ! Fixed parameters of the IMF.
  integer,          parameter                :: imfPieceCount=1
  double precision, dimension(imfPieceCount) :: massLower=[0.1d0],massUpper=[125.0d0],massExponent=[-2.35d0],imfNormalization

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_Salpeter</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_Salpeter(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex=imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_Salpeter

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_Salpeter</unitName>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_Salpeter(imfNames,imfDescriptors)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(inout) :: imfNames(:),imfDescriptors(:)

    imfNames      (imfIndex)="Salpeter"
    imfDescriptors(imfIndex)="Salpeter"
    return
  end subroutine Star_Formation_IMF_Register_Name_Salpeter

  subroutine Star_Formation_IMF_Initialize_Salpeter
    !% Initialize the Salpeter IMF module.
    use Input_Parameters
    use Star_Formation_IMF_PPL
    implicit none

    if (.not.imfSalpeterInitialized) then
       !$omp critical (IMF_Salpeter_Initialize)
       if (.not.imfSalpeterInitialized) then
          !@ <inputParameter>
          !@   <name>imfSalpeterRecycledInstantaneous</name>
          !@   <defaultValue>0.39</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The recycled fraction for the Salpeter \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfSalpeterRecycledInstantaneous',imfSalpeterRecycledInstantaneous,defaultValue=0.39d0)
          !@ <inputParameter>
          !@   <name>imfSalpeterYieldInstantaneous</name>
          !@   <defaultValue>0.02</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The yield for the Salpeter \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfSalpeterYieldInstantaneous'   ,imfSalpeterYieldInstantaneous   ,defaultValue=0.02d0)
          
          ! Get the normalization for this IMF.
          call Piecewise_Power_Law_IMF_Normalize(massLower,massUpper,massExponent,imfNormalization)
          
          imfSalpeterInitialized=.true.
       end if
       !$omp end critical (IMF_Salpeter_Initialize)
    end if
    return
  end subroutine Star_Formation_IMF_Initialize_Salpeter

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_Salpeter</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_Salpeter(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Salpeter
       minimumMass=massLower(1)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_Salpeter

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_Salpeter</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_Salpeter(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Salpeter
       maximumMass=massUpper(imfPieceCount)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_Salpeter

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_Salpeter</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_Salpeter(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(in)    :: initialMass
    double precision, intent(out)   :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Salpeter
       imfPhi=Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_Salpeter

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_Salpeter</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_Salpeter(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Salpeter
       recycledFraction=imfSalpeterRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_Salpeter

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_Salpeter</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_Salpeter(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Salpeter
       yield=imfSalpeterYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_Salpeter

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_Salpeter</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_Salpeter(imfSelected,imfMatched,imfMass,imfPhi)
    !% Register the name of this IMF.
    use Memory_Management
    use Numerical_Ranges
    use Star_Formation_IMF_PPL
    implicit none
    integer,          intent(in)                               :: imfSelected
    logical,          intent(inout)                            :: imfMatched
    double precision, intent(inout), allocatable, dimension(:) :: imfMass,imfPhi
    integer,          parameter                                :: nPoints=100

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Salpeter
       call Alloc_Array(imfMass,[nPoints])
       call Alloc_Array(imfPhi ,[nPoints])
       imfMass=Make_Range(massLower(1),massUpper(imfPieceCount),nPoints,rangeType=rangeTypeLogarithmic)
       imfPhi =Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,imfMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Tabulate_Salpeter

end module Star_Formation_IMF_Salpeter
