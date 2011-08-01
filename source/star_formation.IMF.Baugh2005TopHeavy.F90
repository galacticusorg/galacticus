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


!% Contains a module which implements the top-heavy stellar initial mass function from \cite{baugh_can_2005}.

module Star_Formation_IMF_Baugh2005TopHeavy
  !% Implements the top-heavy stellar initial mass function from \cite{baugh_can_2005}.
  implicit none
  private
  public :: Star_Formation_IMF_Register_Baugh2005TopHeavy, Star_Formation_IMF_Register_Name_Baugh2005TopHeavy,&
       & Star_Formation_IMF_Recycled_Instantaneous_Baugh2005TopHeavy, Star_Formation_IMF_Yield_Instantaneous_Baugh2005TopHeavy,&
       & Star_Formation_IMF_Tabulate_Baugh2005TopHeavy, Star_Formation_IMF_Minimum_Mass_Baugh2005TopHeavy,&
       & Star_Formation_IMF_Maximum_Mass_Baugh2005TopHeavy , Star_Formation_IMF_Phi_Baugh2005TopHeavy

  ! Index assigned to this IMF.
  integer :: imfIndex=-1

  ! Flag indicating if the module has been initialized.
  logical :: imfBaugh2005TopHeavyInitialized=.false.

  ! Parameters of the IMF.
  double precision :: imfBaugh2005TopHeavyRecycledInstantaneous, imfBaugh2005TopHeavyYieldInstantaneous

  ! Fixed parameters of the IMF.
  integer,          parameter                :: imfPieceCount=1
  double precision, dimension(imfPieceCount) :: massLower=[0.15d0],massUpper=[125.0d0],massExponent=[-1.0d0],imfNormalization

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_Baugh2005TopHeavy</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_Baugh2005TopHeavy(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex         =imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_Baugh2005TopHeavy

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_Baugh2005TopHeavy</unitName>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_Baugh2005TopHeavy(imfNames,imfDescriptors)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(inout) :: imfNames(:),imfDescriptors(:)

    imfNames      (imfIndex)="Baugh2005TopHeavy"
    imfDescriptors(imfIndex)="Baugh2005TopHeavy"
    return
  end subroutine Star_Formation_IMF_Register_Name_Baugh2005TopHeavy

  subroutine Star_Formation_IMF_Initialize_Baugh2005TopHeavy
    !% Initialize the Baugh2005TopHeavy IMF module.
    use Input_Parameters
    use Star_Formation_IMF_PPL
    implicit none

    !$omp critical (IMF_Baugh2005TopHeavy_Initialize)
    if (.not.imfBaugh2005TopHeavyInitialized) then
       !@ <inputParameter>
       !@   <name>imfBaugh2005TopHeavyRecycledInstantaneous</name>
       !@   <defaultValue>0.91 \citep{baugh_can_2005}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The recycled fraction for the \cite{baugh_can_2005} top-heavy \IMF\ in the instantaneous recycling approximation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('imfBaugh2005TopHeavyRecycledInstantaneous',imfBaugh2005TopHeavyRecycledInstantaneous,defaultValue=0.57d0)
       !@ <inputParameter>
       !@   <name>imfBaugh2005TopHeavyYieldInstantaneous</name>
       !@   <defaultValue>0.15 \citep{baugh_can_2005}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The yield for the \cite{baugh_can_2005} top-heavy \IMF\ in the instantaneous recycling approximation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('imfBaugh2005TopHeavyYieldInstantaneous'   ,imfBaugh2005TopHeavyYieldInstantaneous   ,defaultValue=0.044d0)

       ! Get the normalization for this IMF.
       call Piecewise_Power_Law_IMF_Normalize(massLower,massUpper,massExponent,imfNormalization)

       imfBaugh2005TopHeavyInitialized=.true.
    end if
    !$omp end critical (IMF_Baugh2005TopHeavy_Initialize)
    return
  end subroutine Star_Formation_IMF_Initialize_Baugh2005TopHeavy

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_Baugh2005TopHeavy</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_Baugh2005TopHeavy(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Baugh2005TopHeavy
       minimumMass=massLower(1)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_Baugh2005TopHeavy

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_Baugh2005TopHeavy</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_Baugh2005TopHeavy(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Baugh2005TopHeavy
       maximumMass=massUpper(imfPieceCount)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_Baugh2005TopHeavy

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_Baugh2005TopHeavy</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_Baugh2005TopHeavy(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(in)    :: initialMass
    double precision, intent(out)   :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Baugh2005TopHeavy
       imfPhi=Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_Baugh2005TopHeavy

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_Baugh2005TopHeavy</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_Baugh2005TopHeavy(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Baugh2005TopHeavy
       recycledFraction=imfBaugh2005TopHeavyRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_Baugh2005TopHeavy

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_Baugh2005TopHeavy</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_Baugh2005TopHeavy(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer,          intent(in)    :: imfSelected
    logical,          intent(inout) :: imfMatched
    double precision, intent(out)   :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Baugh2005TopHeavy
       yield=imfBaugh2005TopHeavyYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_Baugh2005TopHeavy

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_Baugh2005TopHeavy</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_Baugh2005TopHeavy(imfSelected,imfMatched,imfMass,imfPhi)
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
       call Star_Formation_IMF_Initialize_Baugh2005TopHeavy
       call Alloc_Array(imfMass,[nPoints])
       call Alloc_Array(imfPhi ,[nPoints])
       imfMass=Make_Range(massLower(1),massUpper(imfPieceCount),nPoints,rangeType=rangeTypeLogarithmic)
       imfPhi =Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,imfMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Tabulate_Baugh2005TopHeavy

end module Star_Formation_IMF_Baugh2005TopHeavy
