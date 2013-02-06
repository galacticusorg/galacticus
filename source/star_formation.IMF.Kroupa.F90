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

!% Contains a module which implements the Kroupa stellar initial mass function \citep{kroupa_variation_2001}.

module Star_Formation_IMF_Kroupa
  !% Implements the Kroupa stellar initial mass function.
  implicit none
  private
  public :: Star_Formation_IMF_Register_Kroupa, Star_Formation_IMF_Register_Name_Kroupa,&
       & Star_Formation_IMF_Recycled_Instantaneous_Kroupa, Star_Formation_IMF_Yield_Instantaneous_Kroupa,&
       & Star_Formation_IMF_Tabulate_Kroupa, Star_Formation_IMF_Minimum_Mass_Kroupa, Star_Formation_IMF_Maximum_Mass_Kroupa&
       &, Star_Formation_IMF_Phi_Kroupa

  ! Index assigned to this IMF.
  integer :: imfIndex=-1

  ! Flag indicating if the module has been initialized.
  logical :: imfKroupaInitialized=.false.

  ! Parameters of the IMF.
  double precision :: imfKroupaRecycledInstantaneous, imfKroupaYieldInstantaneous

  ! Fixed parameters of the IMF.
  integer,          parameter                :: imfPieceCount=4
  double precision, dimension(imfPieceCount) :: massLower=[0.01d0,0.08d0,0.5d0,1.0d0],massUpper=[0.08d0,0.5d0,1.0d0,125.0d0],massExponent=[&
       &-0.3d0,-1.8d0,-2.7d0,-2.3d0],imfNormalization

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_Kroupa</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_Kroupa(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex=imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_Kroupa

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_Kroupa</unitName>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_Kroupa(imfNames,imfDescriptors)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(inout) :: imfNames(:),imfDescriptors(:)

    imfNames      (imfIndex)="Kroupa"
    imfDescriptors(imfIndex)="Kroupa"
    return
  end subroutine Star_Formation_IMF_Register_Name_Kroupa

  subroutine Star_Formation_IMF_Initialize_Kroupa
    !% Initialize the Kroupa IMF module.
    use Input_Parameters
    use Star_Formation_IMF_PPL
    implicit none

    if (.not.imfKroupaInitialized) then
       !$omp critical (IMF_Kroupa_Initialize)
       if (.not.imfKroupaInitialized) then
          !@ <inputParameter>
          !@   <name>imfKroupaRecycledInstantaneous</name>
          !@   <defaultValue>0.30 (internally computed)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The recycled fraction for the Kroupa \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfKroupaRecycledInstantaneous',imfKroupaRecycledInstantaneous,defaultValue=0.30d0)
          !@ <inputParameter>
          !@   <name>imfKroupaYieldInstantaneous</name>
          !@   <defaultValue>0.023 (internally computed)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The yield for the Kroupa \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfKroupaYieldInstantaneous'   ,imfKroupaYieldInstantaneous   ,defaultValue=0.023d0)
          
          ! Get the normalization for this IMF.
          call Piecewise_Power_Law_IMF_Normalize(massLower,massUpper,massExponent,imfNormalization)
          
          imfKroupaInitialized=.true.
       end if
       !$omp end critical (IMF_Kroupa_Initialize)
    end if
    return
  end subroutine Star_Formation_IMF_Initialize_Kroupa

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_Kroupa</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_Kroupa(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kroupa
       minimumMass=massLower(1)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_Kroupa

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_Kroupa</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_Kroupa(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kroupa
       maximumMass=massUpper(imfPieceCount)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_Kroupa

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_Kroupa</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_Kroupa(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(in)    :: initialMass
    double precision, intent(out)   :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kroupa
       imfPhi=Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_Kroupa

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_Kroupa</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_Kroupa(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kroupa
       recycledFraction=imfKroupaRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_Kroupa

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_Kroupa</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_Kroupa(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kroupa
       yield=imfKroupaYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_Kroupa

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_Kroupa</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_Kroupa(imfSelected,imfMatched,imfMass,imfPhi)
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
       call Star_Formation_IMF_Initialize_Kroupa
       call Alloc_Array(imfMass,[nPoints])
       call Alloc_Array(imfPhi ,[nPoints])
       imfMass=Make_Range(massLower(1),massUpper(imfPieceCount),nPoints,rangeType=rangeTypeLogarithmic)
       imfPhi =Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,imfMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Tabulate_Kroupa

end module Star_Formation_IMF_Kroupa
