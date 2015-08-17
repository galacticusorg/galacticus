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

!% Contains a module which implements the stellar initial mass function used by the \href{http://www.bpass.org.uk/}{BPASS} library.

module Star_Formation_IMF_BPASS
  !% Implements the stellar initial mass function used by the \href{http://www.bpass.org.uk/}{BPASS} library.
  implicit none
  private
  public :: Star_Formation_IMF_Register_BPASS, Star_Formation_IMF_Register_Name_BPASS,&
       & Star_Formation_IMF_Recycled_Instantaneous_BPASS, Star_Formation_IMF_Yield_Instantaneous_BPASS,&
       & Star_Formation_IMF_Tabulate_BPASS, Star_Formation_IMF_Minimum_Mass_BPASS, Star_Formation_IMF_Maximum_Mass_BPASS&
       &, Star_Formation_IMF_Phi_BPASS

  ! Index assigned to this IMF.
  integer                                    :: imfIndex                      =-1

  ! Flag indicating if the module has been initialized.
  logical                                    :: imfBPASSInitialized          =.false.

  ! Parameters of the IMF.
  double precision                           :: imfBPASSRecycledInstantaneous               , imfBPASSYieldInstantaneous

  ! Fixed parameters of the IMF.
  integer         , parameter                :: imfPieceCount                 =2
  double precision, dimension(imfPieceCount) :: imfNormalization                            , massExponent               =[-1.3d0,-2.35d0], &
       &                                        massLower                     =[0.1d0,0.5d0], massUpper                  =[ 0.5d0,120.0d0]

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_BPASS</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_BPASS(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex=imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_BPASS

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_BPASS</unitName>
  !#  <name>BPASS</name>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_BPASS(imfNames,imfDescriptors)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(inout) :: imfDescriptors(:), imfNames(:)

    imfNames      (imfIndex)="BPASS"
    imfDescriptors(imfIndex)="BPASS"
    return
  end subroutine Star_Formation_IMF_Register_Name_BPASS

  subroutine Star_Formation_IMF_Initialize_BPASS
    !% Initialize the BPASS IMF module.
    use Input_Parameters
    use Star_Formation_IMF_PPL
    implicit none

    if (.not.imfBPASSInitialized) then
       !$omp critical (IMF_BPASS_Initialize)
       if (.not.imfBPASSInitialized) then
          !@ <inputParameter>
          !@   <name>imfBPASSRecycledInstantaneous</name>
          !@   <defaultValue>0.30 (internally computed)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The recycled fraction for the BPASS \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfBPASSRecycledInstantaneous',imfBPASSRecycledInstantaneous,defaultValue=0.30d0)
          !@ <inputParameter>
          !@   <name>imfBPASSYieldInstantaneous</name>
          !@   <defaultValue>0.023 (internally computed)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The yield for the BPASS \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfBPASSYieldInstantaneous'   ,imfBPASSYieldInstantaneous   ,defaultValue=0.023d0)

          ! Get the normalization for this IMF.
          call Piecewise_Power_Law_IMF_Normalize(massLower,massUpper,massExponent,imfNormalization)

          imfBPASSInitialized=.true.
       end if
       !$omp end critical (IMF_BPASS_Initialize)
    end if
    return
  end subroutine Star_Formation_IMF_Initialize_BPASS

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_BPASS</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_BPASS(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_BPASS
       minimumMass=massLower(1)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_BPASS

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_BPASS</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_BPASS(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_BPASS
       maximumMass=massUpper(imfPieceCount)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_BPASS

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_BPASS</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_BPASS(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(in   ) :: initialMass
    double precision, intent(  out) :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_BPASS
       imfPhi=Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_BPASS

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_BPASS</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_BPASS(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_BPASS
       recycledFraction=imfBPASSRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_BPASS

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_BPASS</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_BPASS(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_BPASS
       yield=imfBPASSYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_BPASS

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_BPASS</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_BPASS(imfSelected,imfMatched,imf)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    use Tables
    implicit none
    integer                      , intent(in   ) :: imfSelected
    logical                      , intent(inout) :: imfMatched
    class  (table1D), allocatable, intent(inout) :: imf
    integer         , parameter                  :: nPoints    =100

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_BPASS
       allocate(table1DLogarithmicLinear :: imf)
       select type (imf)
       type is (table1DLogarithmicLinear)
          call imf%create  (                                              &
               &            massLower(1            )                    , &
               &            massUpper(imfPieceCount)                    , &
               &            nPoints                                       &
               &           )
          call imf%populate(                                              &
               &            Piecewise_Power_Law_IMF_Phi(                  &
               &                                        massLower       , &
               &                                        massUpper       , &
               &                                        massExponent    , &
               &                                        imfNormalization, &
               &                                        imf%xs()          &
               &                                       )                  &
               &           )
          imfMatched=.true.
       end select
    end if
    return
  end subroutine Star_Formation_IMF_Tabulate_BPASS

end module Star_Formation_IMF_BPASS
