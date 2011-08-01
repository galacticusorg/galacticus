!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements the Kennicutt stellar initial mass function \citep{kennicutt_rate_1983}.

module Star_Formation_IMF_Kennicutt
  !% Implements the Kennicutt stellar initial mass function.
  implicit none
  private
  public :: Star_Formation_IMF_Register_Kennicutt, Star_Formation_IMF_Register_Name_Kennicutt,&
       & Star_Formation_IMF_Recycled_Instantaneous_Kennicutt, Star_Formation_IMF_Yield_Instantaneous_Kennicutt,&
       & Star_Formation_IMF_Tabulate_Kennicutt, Star_Formation_IMF_Minimum_Mass_Kennicutt, Star_Formation_IMF_Maximum_Mass_Kennicutt&
       &, Star_Formation_IMF_Phi_Kennicutt

  ! Index assigned to this IMF.
  integer :: imfIndex=-1

  ! Flag indicating if the module has been initialized.
  logical :: imfKennicuttInitialized=.false.

  ! Parameters of the IMF.
  double precision :: imfKennicuttRecycledInstantaneous, imfKennicuttYieldInstantaneous

  ! Fixed parameters of the IMF.
  integer,          parameter                :: imfPieceCount=3
  double precision, dimension(imfPieceCount) :: massLower=[0.1d0,1.0d0,2.0d0],massUpper=[1.0d0,2.0d0,125.0d0],massExponent=[&
       &-1.25d0,-2.00d0, -2.30d0],imfNormalization

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_Kennicutt</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_Kennicutt(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex=imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_Kennicutt

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_Kennicutt</unitName>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_Kennicutt(imfNames,imfDescriptors)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(inout) :: imfNames(:),imfDescriptors(:)

    imfNames      (imfIndex)="Kennicutt"
    imfDescriptors(imfIndex)="Kennicutt"
    return
  end subroutine Star_Formation_IMF_Register_Name_Kennicutt

  subroutine Star_Formation_IMF_Initialize_Kennicutt
    !% Initialize the Kennicutt IMF module.
    use Input_Parameters
    use Star_Formation_IMF_PPL
    implicit none

    !$omp critical (IMF_Kennicutt_Initialize)
    if (.not.imfKennicuttInitialized) then
       !@ <inputParameter>
       !@   <name>imfKennicuttRecycledInstantaneous</name>
       !@   <defaultValue>0.57 (internally computed)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The recycled fraction for the Kennicutt \IMF\ in the instantaneous recycling approximation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('imfKennicuttRecycledInstantaneous',imfKennicuttRecycledInstantaneous,defaultValue=0.57d0)
       !@ <inputParameter>
       !@   <name>imfKennicuttYieldInstantaneous</name>
       !@   <defaultValue>0.044 (internally computed)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The yield for the Kennicutt \IMF\ in the instantaneous recycling approximation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('imfKennicuttYieldInstantaneous'   ,imfKennicuttYieldInstantaneous   ,defaultValue=0.044d0)

       ! Get the normalization for this IMF.
       call Piecewise_Power_Law_IMF_Normalize(massLower,massUpper,massExponent,imfNormalization)

       imfKennicuttInitialized=.true.
    end if
    !$omp end critical (IMF_Kennicutt_Initialize)
    return
  end subroutine Star_Formation_IMF_Initialize_Kennicutt

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_Kennicutt</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_Kennicutt(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kennicutt
       minimumMass=massLower(1)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_Kennicutt

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_Kennicutt</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_Kennicutt(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kennicutt
       maximumMass=massUpper(imfPieceCount)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_Kennicutt

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_Kennicutt</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_Kennicutt(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(in)    :: initialMass
    double precision, intent(out)   :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kennicutt
       imfPhi=Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_Kennicutt

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_Kennicutt</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_Kennicutt(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kennicutt
       recycledFraction=imfKennicuttRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_Kennicutt

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_Kennicutt</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_Kennicutt(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Kennicutt
       yield=imfKennicuttYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_Kennicutt

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_Kennicutt</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_Kennicutt(imfSelected,imfMatched,imfMass,imfPhi)
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
       call Star_Formation_IMF_Initialize_Kennicutt
       call Alloc_Array(imfMass,[nPoints])
       call Alloc_Array(imfPhi ,[nPoints])
       imfMass=Make_Range(massLower(1),massUpper(imfPieceCount),nPoints,rangeType=rangeTypeLogarithmic)
       imfPhi =Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,imfMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Tabulate_Kennicutt

end module Star_Formation_IMF_Kennicutt
